package stelarx.pro;

import stelarx.Config;
import stelarx.Logging;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Tensor-free S6 scoring. Immutable speciation-child taxon sets are shared;
 * each borrowed state holds only child colour counts and one focused O(n²)
 * insertion slice. No array indexed by all three species is retained.
 */
public final class S6DynamicTripletSupport implements S6ScoreSource {
    private static boolean gpuLoadAttempted, gpuLoaded;
    private static boolean cpuLoadAttempted, cpuLoaded;

    private final int n;
    private int[] parts, partOfChild;
    private int[] children, taxa; // CPU focus fallback only after CUDA upload
    private final int childCount, maxStates;
    private int[] taxonOffsets, taxonChildren;
    private final ArrayBlockingQueue<Session> states;
    private final ConcurrentHashMap<Integer, long[]> guideFocus = new ConcurrentHashMap<>();
    private long gpuIndex, colourIndex;
    private final AtomicLong focusCalls = new AtomicLong(), focusNanos = new AtomicLong();
    private final AtomicLong colourCalls = new AtomicLong(), colourNanos = new AtomicLong();
    private volatile boolean closed;

    private S6DynamicTripletSupport(int n, int[] parts, int[] children, int[] taxa,
                                    int maxStates, long gpuIndex) {
        this.n = n;
        this.parts = parts;
        this.children = children;
        this.taxa = taxa;
        this.childCount = children.length - 1;
        this.maxStates = maxStates;
        this.gpuIndex = gpuIndex;
        partOfChild = new int[children.length - 1];
        for (int p = 0; p + 1 < parts.length; p++)
            for (int child = parts[p]; child < parts[p + 1]; child++)
                partOfChild[child] = p;

        taxonOffsets = new int[n + 1];
        for (int taxon : taxa) taxonOffsets[taxon + 1]++;
        for (int t = 0; t < n; t++) taxonOffsets[t + 1] += taxonOffsets[t];
        taxonChildren = new int[taxa.length];
        int[] cursor = taxonOffsets.clone();
        for (int child = 0; child + 1 < children.length; child++)
            for (int at = children[child]; at < children[child + 1]; at++)
                taxonChildren[cursor[taxa[at]]++] = child;

        states = new ArrayBlockingQueue<>(maxStates);
    }

    public static S6DynamicTripletSupport build(List<Tree> genes, int n, int threads) {
        IntBuilder partOffsets = new IntBuilder();
        IntBuilder childOffsets = new IntBuilder();
        IntBuilder taxonIds = new IntBuilder();
        partOffsets.add(0);
        childOffsets.add(0);
        int[] seen = new int[n];
        int[] epoch = {0};
        for (Tree tree : genes)
            addPartitions(tree.root, tree.postorderArray, seen, epoch,
                partOffsets, childOffsets, taxonIds);
        int[] parts = partOffsets.toArray();
        int[] children = childOffsets.toArray();
        int[] taxa = taxonIds.toArray();
        // Restricted guides can use several workers without keeping one
        // colour state per guide. Full-taxon states remain even more bounded.
        int maxStates = n <= 64
            ? Math.max(2, Math.min(8, Math.min(n, threads)))
            : Math.max(2, Math.min(3, threads));
        S6DynamicTripletSupport support = new S6DynamicTripletSupport(
            n, parts, children, taxa, maxStates, 0);
        boolean wantGpu = Config.getInstance().getComputeMode() == Config.ComputeMode.GPU;
        if (wantGpu && tryLoadGPU()) {
            try {
                // Upload the inverted taxon-to-child index once. Each focused
                // kernel visits only sides containing the inserted species.
                support.gpuIndex = createFocusIndex(n, parts, children, taxa,
                    support.partOfChild, support.taxonOffsets, support.taxonChildren);
                // CUDA owns the immutable child lists now; the Java scorer
                // needs only the inverted membership and partition offsets.
                support.children = null;
                support.taxa = null;
            } catch (IllegalStateException e) {
                if (Config.getInstance().isGpuStrict()) throw e;
                Logging.warn("S6 focused CUDA index unavailable; using CPU: %s", e.getMessage());
            } catch (UnsatisfiedLinkError e) {
                gpuLoaded = false;
                if (Config.getInstance().isGpuStrict())
                    throw new IllegalStateException("S6 focused CUDA entry points unavailable", e);
                Logging.warn("S6 focused CUDA entry points unavailable; using CPU");
            }
        } else if (wantGpu && Config.getInstance().isGpuStrict()) {
            throw new IllegalStateException("S6 focused CUDA library is unavailable");
        }
        if (tryLoadCPU()) {
            try {
                support.colourIndex = createColourIndex(n, parts,
                    support.partOfChild, support.taxonOffsets,
                    support.taxonChildren, maxStates);
            } catch (IllegalStateException | UnsatisfiedLinkError e) {
                Logging.warn("S6 native colour scorer unavailable; using Java: %s",
                    e.getMessage());
            }
        }
        long bytes = 4L * (parts.length
            + (support.children == null ? 0 : children.length + taxa.length)
            + support.partOfChild.length + support.taxonOffsets.length
            + support.taxonChildren.length)
            + (long) maxStates * (3L * support.childCount * Character.BYTES
                + 2L * (parts.length - 1) * Integer.BYTES
                + 2L * n * Integer.BYTES
                + 2L * n * n * Long.BYTES);
        for (int i = 0; i < maxStates; i++)
            support.states.add(support.new Session(i));
        if (support.gpuIndex != 0 && support.colourIndex != 0) {
            // Native CPU and CUDA indexes own their immutable copies. No Java
            // fallback arrays are needed during this GPU search.
            support.parts = support.partOfChild = null;
            support.taxonOffsets = support.taxonChildren = null;
        }
        Logging.info("S6 dynamic support: %d speciation partitions, %d taxa, "
            + "%d MiB indexed/state arrays, %d bounded states (%s focus, %s groups)",
            parts.length - 1, n, bytes >> 20, maxStates,
            support.gpuIndex != 0 ? "CUDA" : "CPU",
            support.colourIndex != 0 ? "native CPU" : "Java CPU");
        return support;
    }

    private static void addPartitions(TreeNode node, int[] leaves,
                                      int[] seen, int[] epoch,
                                      IntBuilder parts, IntBuilder children,
                                      IntBuilder taxa) {
        if (node.isLeaf()) return;
        TreeNode[] sides = node.isPolytomous()
            ? node.children : new TreeNode[]{node.left, node.right};
        for (TreeNode side : sides)
            addPartitions(side, leaves, seen, epoch, parts, children, taxa);
        if (!node.isSpeciation()) return;
        for (TreeNode side : sides) {
            int mark = ++epoch[0];
            int begin = taxa.size();
            for (int pos = side.rangeStart; pos < side.rangeEnd; pos++) {
                int taxon = leaves[pos];
                if (seen[taxon] != mark) {
                    seen[taxon] = mark;
                    taxa.add(taxon);
                }
            }
            taxa.sort(begin);
            children.add(taxa.size());
        }
        parts.add(children.size() - 1);
    }

    @Override public int taxa() { return n; }
    public int scoringStates() { return maxStates; }

    @Override public S6ScoreSource fork() {
        if (closed) throw new IllegalStateException("S6 support is closed");
        try {
            Session state = states.take();
            state.borrowed = true;
            return state;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("Interrupted waiting for S6 scorer", e);
        }
    }

    @Override public Focus focus(int taxon) {
        throw new IllegalStateException("Borrow an S6 scoring state first");
    }

    @Override public void prepareGuideFocus(int[] selected) {
        // All guides in one round insert the same species in different orders.
        // Reuse only when the entire batch fits a fixed 32-MiB cache; large-n
        // guide sets therefore cannot recreate cubic support storage.
        if ((long) selected.length * 2L * n * n * Long.BYTES > 32L << 20)
            return;
        for (int taxon : selected) {
            if (guideFocus.containsKey(taxon)) continue;
            long[] slice = new long[2 * n * n];
            if (gpuIndex != 0) fillFocus(gpuIndex, taxon, slice);
            else fillFocusCPU(taxon, slice);
            guideFocus.put(taxon, slice);
        }
    }

    @Override public void clearGuideFocus() { guideFocus.clear(); }
    @Override public long[] groupScores(int[] a, int[] b, int[] c) {
        throw new IllegalStateException("Borrow an S6 scoring state first");
    }
    @Override public long splitWeight(int[] a, int[] b) {
        throw new IllegalStateException("Borrow an S6 scoring state first");
    }

    @Override public synchronized void close() {
        if (closed) return;
        closed = true;
        guideFocus.clear();
        Logging.info("S6 dynamic scorer: %d focused builds in %d ms, "
            + "%d colour queries in %d ms", focusCalls.get(),
            focusNanos.get() / 1_000_000, colourCalls.get(),
            colourNanos.get() / 1_000_000);
        if (gpuIndex != 0) destroyFocusIndex(gpuIndex);
        if (colourIndex != 0) destroyColourIndex(colourIndex);
    }

    private static synchronized boolean tryLoadGPU() {
        if (!gpuLoadAttempted) {
            gpuLoadAttempted = true;
            try {
                System.loadLibrary("stelar_pro_s6_triplets_gpu");
                gpuLoaded = true;
            } catch (UnsatisfiedLinkError ignored) {
                gpuLoaded = false;
            }
        }
        return gpuLoaded;
    }

    private static synchronized boolean tryLoadCPU() {
        if (!cpuLoadAttempted) {
            cpuLoadAttempted = true;
            try {
                System.loadLibrary("stelar_pro_s6_triplets");
                cpuLoaded = true;
            } catch (UnsatisfiedLinkError ignored) {
                cpuLoaded = false;
            }
        }
        return cpuLoaded;
    }

    private static native long createFocusIndex(int n, int[] parts,
                                                int[] children, int[] taxa,
                                                int[] partOfChild, int[] taxonOffsets,
                                                int[] taxonChildren);
    private static native void fillFocus(long index, int taxon, long[] destination);
    private static native void destroyFocusIndex(long index);
    private static native long createColourIndex(int n, int[] parts,
                                                int[] partOfChild,
                                                int[] taxonOffsets,
                                                int[] taxonChildren,
                                                int stateCount);
    private static native void scoreGroupsNative(long index, int stateId,
                                                 int[] desired, long[] scores);
    private static native long scoreSplitNative(long index, int stateId,
                                                int[] desired);
    private static native void destroyColourIndex(long index);

    private void fillFocusCPU(int x, long[] result) {
        Arrays.fill(result, 0L);
        int paired = n * n;
        for (int p = 0; p + 1 < parts.length; p++) {
            for (int side = parts[p]; side < parts[p + 1]; side++) {
                if (Arrays.binarySearch(taxa, children[side], children[side + 1], x) < 0)
                    continue;
                for (int out = parts[p]; out < parts[p + 1]; out++) {
                    if (out == side) continue;
                    for (int i = children[out]; i < children[out + 1]; i++)
                        for (int j = i + 1; j < children[out + 1]; j++) {
                            int a = taxa[i], b = taxa[j];
                            if (a != x && b != x) result[a * n + b]++;
                        }
                    for (int i = children[out]; i < children[out + 1]; i++) {
                        int c = taxa[i];
                        if (c == x) continue;
                        for (int j = children[side]; j < children[side + 1]; j++) {
                            int a = taxa[j];
                            if (a != x && a != c) result[paired + c * n + a]++;
                        }
                    }
                }
            }
        }
    }

    /** One mutable colour state; states are reused rather than copied per guide. */
    private final class Session implements S6ScoreSource {
        // n is bounded by the focused-array index limit, so a distinct-species
        // count always fits in an unsigned 16-bit Java char.
        private final char[] colourCounts;
        private final int[] colours, desired = new int[n];
        private final int[] partMarks, touched;
        private final long[] focusData = new long[2 * n * n];
        private final long[] totals;
        private final int stateId;
        private int epoch, lastFocusTaxon = -1;
        private boolean borrowed;

        private Session(int stateId) {
            this.stateId = stateId;
            if (colourIndex == 0) {
                colourCounts = new char[3 * childCount];
                colours = new int[n];
                Arrays.fill(colours, -1);
                partMarks = new int[parts.length - 1];
                touched = new int[parts.length - 1];
                totals = new long[3];
            } else {
                colourCounts = null;
                colours = partMarks = touched = null;
                totals = null;
            }
        }

        @Override public int taxa() { return n; }
        @Override public S6ScoreSource fork() { return S6DynamicTripletSupport.this.fork(); }

        @Override public Focus focus(int x) {
            long[] slice = guideFocus.get(x);
            if (slice == null) {
                slice = focusData;
                if (lastFocusTaxon != x) {
                    long started = System.nanoTime();
                    if (gpuIndex != 0) fillFocus(gpuIndex, x, slice);
                    else fillFocusCPU(x, slice);
                    lastFocusTaxon = x;
                    focusCalls.incrementAndGet();
                    focusNanos.addAndGet(System.nanoTime() - started);
                }
            }
            final long[] data = slice;
            return (a, b, c) -> {
                if (a == b || a == c || b == c) return 0;
                if (c == x) {
                    if (a > b) { int swap = a; a = b; b = swap; }
                    return data[a * n + b];
                }
                int other = a == x ? b : b == x ? a : -1;
                if (other < 0)
                    throw new IllegalArgumentException("Triplet does not contain focused taxon");
                return data[n * n + c * n + other];
            };
        }

        @Override public long[] groupScores(int[] a, int[] b, int[] c) {
            long started = System.nanoTime();
            Arrays.fill(desired, -1);
            for (int taxon : a) desired[taxon] = 0;
            for (int taxon : b) desired[taxon] = 1;
            for (int taxon : c) desired[taxon] = 2;
            long[] result;
            if (colourIndex != 0) {
                result = new long[3];
                scoreGroupsNative(colourIndex, stateId, desired, result);
            } else {
                updateColours();
                result = totals.clone();
            }
            colourCalls.incrementAndGet();
            colourNanos.addAndGet(System.nanoTime() - started);
            return result;
        }

        @Override public long splitWeight(int[] a, int[] b) {
            Arrays.fill(desired, -1);
            for (int taxon : a) desired[taxon] = 0;
            for (int taxon : b) desired[taxon] = 1;
            if (colourIndex != 0)
                return scoreSplitNative(colourIndex, stateId, desired);
            updateColours();
            long score = 0;
            for (int p = 0; p + 1 < parts.length; p++) {
                int sumA = 0, sumB = 0;
                for (int child = parts[p]; child < parts[p + 1]; child++) {
                    sumA += colourCounts[3 * child];
                    sumB += colourCounts[3 * child + 1];
                }
                for (int child = parts[p]; child < parts[p + 1]; child++) {
                    long ca = colourCounts[3 * child];
                    long cb = colourCounts[3 * child + 1];
                    score += ca * (ca - 1) / 2 * (sumB - cb)
                           + cb * (cb - 1) / 2 * (sumA - ca);
                }
            }
            return score;
        }

        /** Recolour only gene partitions containing changed species. */
        private void updateColours() {
            if (++epoch == Integer.MAX_VALUE) {
                Arrays.fill(partMarks, 0);
                epoch = 1;
            }
            int touchedCount = 0;
            for (int taxon = 0; taxon < n; taxon++) {
                if (desired[taxon] == colours[taxon]) continue;
                for (int at = taxonOffsets[taxon]; at < taxonOffsets[taxon + 1]; at++) {
                    int part = partOfChild[taxonChildren[at]];
                    if (partMarks[part] != epoch) {
                        partMarks[part] = epoch;
                        touched[touchedCount++] = part;
                    }
                }
            }
            for (int i = 0; i < touchedCount; i++) accumulate(touched[i], -1);
            for (int taxon = 0; taxon < n; taxon++) {
                int before = colours[taxon], after = desired[taxon];
                if (before == after) continue;
                for (int at = taxonOffsets[taxon]; at < taxonOffsets[taxon + 1]; at++) {
                    int base = 3 * taxonChildren[at];
                    if (before >= 0) colourCounts[base + before]--;
                    if (after >= 0) colourCounts[base + after]++;
                }
                colours[taxon] = after;
            }
            for (int i = 0; i < touchedCount; i++) accumulate(touched[i], 1);
        }

        /** Each side pair is counted once, matching the old tensor builder. */
        private void accumulate(int p, int sign) {
            long sumA = 0, sumB = 0, sumC = 0;
            for (int child = parts[p]; child < parts[p + 1]; child++) {
                int base = 3 * child;
                sumA += colourCounts[base];
                sumB += colourCounts[base + 1];
                sumC += colourCounts[base + 2];
            }
            for (int child = parts[p]; child < parts[p + 1]; child++) {
                int base = 3 * child;
                long a = colourCounts[base], b = colourCounts[base + 1];
                long c = colourCounts[base + 2];
                totals[0] += sign * a * b * (sumC - c);
                totals[1] += sign * a * c * (sumB - b);
                totals[2] += sign * b * c * (sumA - a);
            }
        }

        @Override public void close() {
            if (!borrowed) return;
            borrowed = false;
            states.add(this);
        }
    }

    private static final class IntBuilder {
        private int[] values = new int[1024];
        private int size;
        int size() { return size; }
        void add(int value) {
            if (size == values.length) values = Arrays.copyOf(values, values.length * 2);
            values[size++] = value;
        }
        void sort(int begin) { Arrays.sort(values, begin, size); }
        int[] toArray() { return Arrays.copyOf(values, size); }
    }
}
