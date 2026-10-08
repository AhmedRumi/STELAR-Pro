package stelarx.pro;

import stelarx.Logging;
import stelarx.Config;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.Arrays;
import java.util.List;

/**
 * Exact rooted-triplet support of the original tagged gene families.
 *
 * <p>{@code support(a,b,c)} counts the speciation-rooted gene partitions
 * supporting {@code ab|c}, with each species represented once in each child.
 * Consequently, summing these counts at a species-tree split reproduces the
 * same binary-partition formula as {@code WeightTable}. No guide tree is an
 * observation, and duplication-rooted gene nodes never contribute.</p>
 */
public final class S6TripletSupport {
    private static boolean loadAttempted;
    private static boolean loaded;
    private static boolean gpuLoadAttempted;
    private static boolean gpuLoaded;

    private final int n;
    /** Outgroup-major; only entries with a < b and a,b != outgroup are used. */
    private final long[] counts;
    private final int partitionCount;

    private S6TripletSupport(int n, long[] counts, int partitionCount) {
        this.n = n;
        this.counts = counts;
        this.partitionCount = partitionCount;
    }

    public int taxa() { return n; }
    public int partitionCount() { return partitionCount; }

    public long support(int a, int b, int outgroup) {
        if (a == b || a == outgroup || b == outgroup) return 0L;
        if (a > b) { int t = a; a = b; b = t; }
        return counts[(outgroup * n + a) * n + b];
    }

    /** Score one species split; useful as an independent parity oracle. */
    public long splitWeight(int[] a, int[] b) {
        long score = 0;
        for (int i = 0; i < a.length; i++)
            for (int j = i + 1; j < a.length; j++)
                for (int taxon : b) score += support(a[i], a[j], taxon);
        for (int i = 0; i < b.length; i++)
            for (int j = i + 1; j < b.length; j++)
                for (int taxon : a) score += support(b[i], b[j], taxon);
        return score;
    }

    public static S6TripletSupport build(List<Tree> genes, int n, int threads) {
        if (n < 3 || (long) n * n * n > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("S6 triplet tensor exceeds Java array limits");
        }
        IntBuilder partOffsets = new IntBuilder();
        IntBuilder childOffsets = new IntBuilder();
        IntBuilder taxa = new IntBuilder();
        partOffsets.add(0);
        childOffsets.add(0);
        int[] seen = new int[n];
        int[] epoch = {0};
        for (Tree gene : genes) {
            addPartitions(gene.root, gene.postorderArray, seen, epoch,
                partOffsets, childOffsets, taxa);
        }
        int parts = partOffsets.size() - 1;
        long started = System.nanoTime();
        long[] counts;
        String backend;
        boolean requestedGPU = Config.getInstance().getComputeMode()
            == Config.ComputeMode.GPU;
        if (requestedGPU && !tryLoadGPU()
                && Config.getInstance().isGpuStrict()) {
            throw new IllegalStateException(
                "S6 CUDA triplet-support library is unavailable under --gpu-strict");
        }
        if (requestedGPU && gpuLoaded) {
            try {
                counts = buildNativeGPU(n, partOffsets.toArray(),
                    childOffsets.toArray(), taxa.toArray());
                backend = "CUDA";
            } catch (IllegalStateException e) {
                if (Config.getInstance().isGpuStrict()) throw e;
                Logging.warn("S6 CUDA support construction failed; using CPU: %s",
                    e.getMessage());
                counts = buildCPU(n, partOffsets, childOffsets, taxa, threads);
                backend = loaded ? "native CPU fallback" : "Java CPU fallback";
            }
        } else {
            if (requestedGPU) Logging.warn(
                "S6 CUDA triplet-support library unavailable; using CPU construction");
            counts = buildCPU(n, partOffsets, childOffsets, taxa, threads);
            backend = loaded ? "native CPU" : "Java CPU";
        }
        Logging.info("S6 triplet support: %d speciation partitions, %d taxa, "
            + "%d MiB, %d ms (%s)", parts, n, counts.length * Long.BYTES / (1 << 20),
            (System.nanoTime() - started) / 1_000_000, backend);
        return new S6TripletSupport(n, counts, parts);
    }

    private static long[] buildCPU(int n, IntBuilder partOffsets,
                                   IntBuilder childOffsets, IntBuilder taxa,
                                   int threads) {
        if (tryLoad()) {
            return buildNative(n, partOffsets.toArray(), childOffsets.toArray(),
                taxa.toArray(), Math.max(1, threads));
        }
        // Small tests and development environments can still exercise S6.
        return buildJava(n, partOffsets.toArray(), childOffsets.toArray(),
            taxa.toArray());
    }

    private static void addPartitions(TreeNode node, int[] leaves,
                                      int[] seen, int[] epoch,
                                      IntBuilder partOffsets,
                                      IntBuilder childOffsets,
                                      IntBuilder taxa) {
        if (node.isLeaf()) return;
        TreeNode[] children = node.isPolytomous()
            ? node.children : new TreeNode[]{node.left, node.right};
        for (TreeNode child : children)
            addPartitions(child, leaves, seen, epoch, partOffsets, childOffsets, taxa);
        if (!node.isSpeciation()) return;

        for (TreeNode child : children) {
            int mark = ++epoch[0];
            int begin = taxa.size();
            for (int pos = child.rangeStart; pos < child.rangeEnd; pos++) {
                int taxon = leaves[pos];
                if (seen[taxon] != mark) {
                    seen[taxon] = mark;
                    taxa.add(taxon);
                }
            }
            taxa.sort(begin);
            childOffsets.add(taxa.size());
        }
        partOffsets.add(childOffsets.size() - 1);
    }

    private static synchronized boolean tryLoad() {
        if (!loadAttempted) {
            loadAttempted = true;
            try {
                System.loadLibrary("stelar_pro_s6_triplets");
                loaded = true;
            } catch (UnsatisfiedLinkError e) {
                loaded = false;
            }
        }
        return loaded;
    }

    private static synchronized boolean tryLoadGPU() {
        if (!gpuLoadAttempted) {
            gpuLoadAttempted = true;
            try {
                System.loadLibrary("stelar_pro_s6_triplets_gpu");
                gpuLoaded = true;
            } catch (UnsatisfiedLinkError e) {
                gpuLoaded = false;
            }
        }
        return gpuLoaded;
    }

    private static native long[] buildNative(int n, int[] partOffsets,
                                             int[] childOffsets, int[] taxa,
                                             int threads);
    private static native long[] buildNativeGPU(int n, int[] partOffsets,
                                                int[] childOffsets, int[] taxa);

    private static long[] buildJava(int n, int[] partOffsets,
                                    int[] childOffsets, int[] taxa) {
        long[] tensor = new long[n * n * n];
        for (int c = 0; c < n; c++) {
            for (int p = 0; p + 1 < partOffsets.length; p++) {
                for (int out = partOffsets[p]; out < partOffsets[p + 1]; out++) {
                    if (Arrays.binarySearch(taxa, childOffsets[out],
                            childOffsets[out + 1], c) < 0) continue;
                    for (int side = partOffsets[p]; side < partOffsets[p + 1]; side++) {
                        if (side == out) continue;
                        for (int i = childOffsets[side]; i < childOffsets[side + 1]; i++) {
                            int a = taxa[i];
                            if (a == c) continue;
                            for (int j = i + 1; j < childOffsets[side + 1]; j++) {
                                int b = taxa[j];
                                if (b != c) tensor[(c * n + a) * n + b]++;
                            }
                        }
                    }
                }
            }
        }
        return tensor;
    }

    private static final class IntBuilder {
        private int[] values = new int[1024];
        private int size;

        int size() { return size; }
        void add(int value) {
            if (size == values.length) values = Arrays.copyOf(values, values.length * 2);
            values[size++] = value;
        }
        void sort(int from) { Arrays.sort(values, from, size); }
        int[] toArray() { return Arrays.copyOf(values, size); }
    }
}
