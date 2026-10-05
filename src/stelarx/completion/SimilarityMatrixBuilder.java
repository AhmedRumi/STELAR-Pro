package stelarx.completion;

import stelarx.Config;
import stelarx.Logging;
import stelarx.gpu.GPUSimilarityMatrix;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.ProgressBar;
import stelarx.util.Threading;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Builds the quartet-based taxon similarity matrix from gene trees.
 *
 * Reproduces ASTRAL-MP's SimilarityMatrix.populateByQuartetDistance exactly:
 * for each tree T with kt leaves, and for each pair (a, b) both present in T,
 *   num_T(a, b) = same_side_T(a, b)
 *              = number of quartets {a, b, x, y} resolved by T where a and b
 *                fall on the same side of the bipartition
 *   den_T(a, b) = C2(kt − 2)
 * The final similarity is sim(a,b) = Σ num_T / Σ den_T.
 *
 * ── CPU path (buildCPU) ──────────────────────────────────────────────────────
 * Per-tree scatter, identical to ASTRAL-MP:
 *   For every internal node u, form three "components":
 *     left  = leaves of left  child of u
 *     right = leaves of right child of u
 *     others = leaves in T but not in subtree(u)
 *   totalPairs = C2(|left|) + C2(|right|) + C2(|others|)
 *   For each component-pair (X, Y), the per-(l ∈ X, r ∈ Y) scatter is
 *     sim = totalPairs − C2(|X|) − C2(|Y|)
 *
 * This sums contributions across every internal node on the unrooted path
 * between a and b — exactly the same-side quartet count. The earlier
 * implementation only scattered at the LCA, dropping the (path × others)
 * contributions; that was the source of the systematic mismatch.
 *
 * ── GPU path (buildGPU) ──────────────────────────────────────────────────────
 * Uses the validated bridge identity
 *   same_side_T(a, b)  =  C2(kt − 2)  −  QD_gt(a, b)
 * with an O(1) closed form for QD_gt via Euler tour + sparse-table RMQ
 * carrying (s, F) child-of-LCA payloads. See
 *   DOCS/similarity-matrix-design.md
 * for the derivation. Ordinary tours use the compact unsigned-16 RMQ; larger
 * tours use an exact blocked RMQ with wide positions and child sizes. The CPU
 * path remains the byte-for-byte-correct reference.
 */
public class SimilarityMatrixBuilder {

    private static final int COMPACT_MAX_TOUR = Character.MAX_VALUE + 1;
    private static final int MAX_JAVA_ARRAY_LENGTH = Integer.MAX_VALUE - 8;
    /**
     * Host flattening cap used only when the established one-shot wide layout
     * cannot fit. Native CUDA batching is normally 512 MiB, so 1 GiB keeps at
     * most two ordinary GPU batches resident on the Java side without creating
     * excessive JNI/setup overhead.
     */
    private static final long STREAMED_WIDE_HOST_BATCH_BYTES = 1L << 30;
    /** Large-N only: enough to keep the observed 50k×1000 wide tree data in one GPU batch. */
    private static final int LARGE_PACKED_GPU_TREE_CAP_MIB = 8192;

    // ── Public entry points ───────────────────────────────────────────────────

    public static SimilarityMatrix buildCPU(List<Tree> trees, int n) {
        preflightPackedHeap(n);
        SimilarityMatrix sm = new SimilarityMatrix(n);

        ProgressBar bar = new ProgressBar("Building similarity matrix (CPU)", trees.size());
        int done = 0;
        for (Tree tree : trees) {
            accumulateCPU(tree, sm);
            bar.update(++done);
        }
        bar.done();

        sm.normalize();
        return sm;
    }

    /**
     * Build S4's rooted-triplet sister-frequency matrix on the CPU.
     *
     * <p>For a pair present in a single-copy rooted tree with {@code k} leaves,
     * {@code k - size(LCA(a,b))} third taxa display {@code ab|c}, out of the
     * {@code k - 2} possible third taxa. Rows are pair-owned, so CPU workers
     * write disjoint matrix cells and need no locks or thread-local matrices.
     */
    public static SimilarityMatrix buildTripletCPU(List<Tree> trees, int n) {
        preflightPackedHeap(n);
        SimilarityMatrix sm = new SimilarityMatrix(n);
        int eMax = maximumTourLength(trees);
        boolean wide = eMax > COMPACT_MAX_TOUR;

        Logging.info("Building rooted-triplet tour data for %d trees (parallel CPU, %s RMQ)",
            trees.size(), wide ? "wide blocked" : "compact");
        ProgressBar tourBar = new ProgressBar("Triplet tour + RMQ build", trees.size());
        AtomicInteger toursDone = new AtomicInteger();

        if (!wide) {
            TripletTourBuilder.CompactData[] tours =
                new TripletTourBuilder.CompactData[trees.size()];
            Threading.processRangeParallel(trees.size(), i -> {
                try {
                    tours[i] = TripletTourBuilder.buildCompact(trees.get(i), n);
                    tourBar.update(toursDone.incrementAndGet());
                } catch (RuntimeException e) {
                    throw treeBuildFailure(i, trees.get(i), e);
                }
            });
            tourBar.done();
            accumulateTripletPairsCompact(tours, sm);
        } else {
            TripletTourBuilder.WideData[] tours =
                new TripletTourBuilder.WideData[trees.size()];
            Threading.processRangeParallel(trees.size(), i -> {
                try {
                    tours[i] = TripletTourBuilder.buildWide(trees.get(i), n);
                    tourBar.update(toursDone.incrementAndGet());
                } catch (RuntimeException e) {
                    throw treeBuildFailure(i, trees.get(i), e);
                }
            });
            tourBar.done();
            accumulateTripletPairsWide(tours, sm);
        }
        sm.normalize();
        return sm;
    }

    /** GPU equivalent of {@link #buildTripletCPU}; uses the existing tiling and batching policy. */
    public static SimilarityMatrix buildTripletGPU(List<Tree> trees, int n) {
        preflightPackedHeap(n);
        int eMax = maximumTourLength(trees);
        int logMax = sparseLog(eMax);
        boolean compactOverflow = eMax > COMPACT_MAX_TOUR;
        boolean compactLayoutOverflow = !fitsJavaArray((long)trees.size() * nextPowerOfTwo(eMax))
            || !fitsJavaArray((long)trees.size() * logMax * nextPowerOfTwo(eMax))
            || !fitsJavaArray((long)trees.size() * n);
        boolean forceWide = Boolean.getBoolean("stelarpro.similarity.forceWide");
        if (forceWide || compactOverflow || compactLayoutOverflow) {
            Logging.info("GPU rooted-triplet similarity RMQ: wide blocked mode");
            return buildTripletGPUWide(trees, n, eMax);
        }
        Logging.debug("GPU rooted-triplet similarity RMQ: compact unsigned-16 mode");
        return buildTripletGPUCompact(trees, n);
    }

    private static void accumulateTripletPairsCompact(
            TripletTourBuilder.CompactData[] tours, SimilarityMatrix sm) {
        ProgressBar bar = new ProgressBar("Building triplet similarity matrix (CPU)", sm.n);
        AtomicInteger done = new AtomicInteger();
        Threading.processRangeParallel(sm.n, a -> {
            for (int b = a + 1; b < sm.n; b++) {
                double numerator = 0.0;
                double denominator = 0.0;
                for (TripletTourBuilder.CompactData tour : tours) {
                    int fa = tour.firstOcc[a], fb = tour.firstOcc[b];
                    if (fa < 0 || fb < 0 || tour.leafCount < 3) continue;
                    int lo = Math.min(fa, fb), hi = Math.max(fa, fb);
                    int lca = TripletTourBuilder.queryCompact(tour, lo, hi);
                    numerator += tour.leafCount - Short.toUnsignedInt(tour.nodeSizes[lca]);
                    denominator += tour.leafCount - 2;
                }
                storePair(sm, a, b, numerator, denominator);
            }
            bar.update(done.incrementAndGet());
        });
        bar.done();
    }

    private static void accumulateTripletPairsWide(
            TripletTourBuilder.WideData[] tours, SimilarityMatrix sm) {
        ProgressBar bar = new ProgressBar("Building triplet similarity matrix (CPU)", sm.n);
        AtomicInteger done = new AtomicInteger();
        Threading.processRangeParallel(sm.n, a -> {
            for (int b = a + 1; b < sm.n; b++) {
                double numerator = 0.0;
                double denominator = 0.0;
                for (TripletTourBuilder.WideData tour : tours) {
                    int fa = tour.firstOcc[a], fb = tour.firstOcc[b];
                    if (fa < 0 || fb < 0 || tour.leafCount < 3) continue;
                    int lo = Math.min(fa, fb), hi = Math.max(fa, fb);
                    int lca = TripletTourBuilder.queryWide(tour, lo, hi);
                    numerator += tour.leafCount - tour.nodeSizes[lca];
                    denominator += tour.leafCount - 2;
                }
                storePair(sm, a, b, numerator, denominator);
            }
            bar.update(done.incrementAndGet());
        });
        bar.done();
    }

    private static void storePair(SimilarityMatrix sm, int a, int b,
                                  double numerator, double denominator) {
        if (denominator == 0.0) return;
        if (sm.isPacked()) {
            sm.addPackedNumerator(a, b, numerator);
            sm.addPackedDenominator(a, b, denominator);
        } else {
            int ab = a * sm.n + b;
            int ba = b * sm.n + a;
            sm.numSum[ab] = sm.numSum[ba] = numerator;
            sm.denSum[ab] = sm.denSum[ba] = denominator;
        }
    }

    private static SimilarityMatrix buildTripletGPUCompact(List<Tree> trees, int n) {
        int k = trees.size();
        TripletTourBuilder.CompactData[] tours = new TripletTourBuilder.CompactData[k];
        ProgressBar buildBar = new ProgressBar("Triplet tour + RMQ build", k);
        AtomicInteger built = new AtomicInteger();
        Threading.processRangeParallel(k, i -> {
            try {
                tours[i] = TripletTourBuilder.buildCompact(trees.get(i), n);
                buildBar.update(built.incrementAndGet());
            } catch (RuntimeException e) {
                throw treeBuildFailure(i, trees.get(i), e);
            }
        });
        buildBar.done();

        int eMaxRaw = 1, logMax = 1;
        for (TripletTourBuilder.CompactData tour : tours) {
            eMaxRaw = Math.max(eMaxRaw, tour.tourLen);
            logMax = Math.max(logMax, tour.log);
        }
        int eMax = nextPowerOfTwo(eMaxRaw);
        int eulerSize = checkedLength((long)k * eMax, "triplet Euler arrays");
        int sparseSize = checkedLength((long)k * logMax * eMax, "triplet sparse RMQ");
        int leafSize = checkedLength((long)k * n, "triplet first-occurrence array");
        Logging.info("Building compact rooted-triplet GPU data: E_max=%d LOG=%d (%.1f MB)",
            eMax, logMax,
            (eulerSize * 2.0 * Short.BYTES + sparseSize * Character.BYTES
                + leafSize * Integer.BYTES) / 1e6);

        short[] depths = new short[eulerSize];
        short[] nodeSizes = new short[eulerSize];
        char[] argmin = new char[sparseSize];
        int[] firstOcc = new int[leafSize];
        int[] eulerLen = new int[k];
        int[] leafCount = new int[k];
        Arrays.fill(firstOcc, -1);
        ProgressBar flatBar = new ProgressBar("Flattening triplet tour data", k);
        AtomicInteger flattened = new AtomicInteger();
        final int flatLogMax = logMax;
        Threading.processRangeParallel(k, i -> {
            TripletTourBuilder.CompactData tour = tours[i];
            int eulerOff = i * eMax;
            System.arraycopy(tour.depths, 0, depths, eulerOff, tour.tourLen);
            System.arraycopy(tour.nodeSizes, 0, nodeSizes, eulerOff, tour.tourLen);
            int sparseOff = i * flatLogMax * eMax;
            for (int level = 0; level < tour.log; level++) {
                int rowLen = Math.max(0, tour.tourLen - (1 << level) + 1);
                System.arraycopy(tour.sparseArgmin[level], 0, argmin,
                    sparseOff + level * eMax, rowLen);
            }
            System.arraycopy(tour.firstOcc, 0, firstOcc, i * n, n);
            eulerLen[i] = tour.tourLen;
            leafCount[i] = tour.leafCount;
            flatBar.update(flattened.incrementAndGet());
        });
        flatBar.done();

        Config cfg = Config.getInstance();
        SimilarityMatrix sm = new SimilarityMatrix(n);
        GPUSimilarityMatrix.computeSimilarityGPU(
            depths, nodeSizes, true,
            null, null, null, null, null,
            argmin, firstOcc, eulerLen, leafCount,
            k, n, eMax, logMax,
            cfg.getGpuDistTileSizeB(), effectiveTreeCapMiB(sm, cfg),
            cfg.getGpuDpProgressInterval(), cfg.getGpuDpProgressMaxSteps(),
            sm.numSum, sm.denSum,
            sm.isPacked() ? sm.packedNumeratorSegments() : null,
            sm.isPacked() ? sm.packedDenominatorSegments() : null,
            sm.packedSegmentShift());
        sm.normalize();
        return sm;
    }

    private static SimilarityMatrix buildTripletGPUWide(List<Tree> trees, int n, int eMaxRaw) {
        if (!wideLayoutFits(trees.size(), n, eMaxRaw, MAX_JAVA_ARRAY_LENGTH)) {
            int batchTrees = wideBatchTreeCount(trees.size(), n, eMaxRaw,
                MAX_JAVA_ARRAY_LENGTH, STREAMED_WIDE_HOST_BATCH_BYTES);
            SimilarityMatrix sm = new SimilarityMatrix(n);
            sm.markStreamedHostBatches();
            int batches = (trees.size() + batchTrees - 1) / batchTrees;
            for (int from = 0, batch = 1; from < trees.size(); batch++) {
                int to = Math.min(trees.size(), from + batchTrees);
                Logging.info("Wide triplet-similarity host batch %d/%d: trees %d..%d",
                    batch, batches, from, to - 1);
                accumulateTripletGPUWide(trees.subList(from, to), n, eMaxRaw, sm, from);
                from = to;
            }
            sm.normalize();
            return sm;
        }
        SimilarityMatrix sm = accumulateTripletGPUWide(trees, n, eMaxRaw, null, 0);
        sm.normalize();
        return sm;
    }

    private static SimilarityMatrix accumulateTripletGPUWide(
            List<Tree> trees, int n, int eMaxRaw, SimilarityMatrix sm, int treeIndexBase) {
        int k = trees.size();
        int eMax = Math.max(1, eMaxRaw);
        int blockSize = EulerTourBuilder.WIDE_BLOCK_SIZE;
        int microLog = EulerTourBuilder.WIDE_MICRO_LOG;
        int blockMax = (eMax + blockSize - 1) / blockSize;
        int macroLog = sparseLog(blockMax);
        int eulerSize = checkedLength((long)k * eMax, "wide triplet Euler arrays");
        int microSize = checkedLength((long)k * microLog * eMax, "wide triplet micro RMQ");
        int macroSize = checkedLength((long)k * macroLog * blockMax, "wide triplet macro RMQ");
        int leafSize = checkedLength((long)k * n, "wide triplet first-occurrence array");

        int[] depths = new int[eulerSize];
        int[] nodeSizes = new int[eulerSize];
        byte[] micro = new byte[microSize];
        int[] macro = new int[macroSize];
        int[] firstOcc = new int[leafSize];
        int[] eulerLen = new int[k];
        int[] leafCount = new int[k];
        Arrays.fill(firstOcc, -1);
        Logging.info("Building wide rooted-triplet GPU data for %d trees (parallel CPU)", k);
        ProgressBar bar = new ProgressBar("Wide triplet tour + blocked RMQ build", k);
        AtomicInteger done = new AtomicInteger();
        Threading.processRangeParallel(k, i -> {
            try {
                TripletTourBuilder.WideData tour =
                    TripletTourBuilder.buildWide(trees.get(i), n);
                int eulerOff = i * eMax;
                System.arraycopy(tour.depths, 0, depths, eulerOff, tour.tourLen);
                System.arraycopy(tour.nodeSizes, 0, nodeSizes, eulerOff, tour.tourLen);
                int microOff = i * microLog * eMax;
                for (int level = 0; level < microLog; level++) {
                    System.arraycopy(tour.microArgmin[level], 0, micro,
                        microOff + level * eMax, tour.tourLen);
                }
                int macroOff = i * macroLog * blockMax;
                for (int level = 0; level < tour.macroLog; level++) {
                    int rowLen = Math.max(0, tour.blockCount - (1 << level) + 1);
                    System.arraycopy(tour.macroArgmin[level], 0, macro,
                        macroOff + level * blockMax, rowLen);
                }
                System.arraycopy(tour.firstOcc, 0, firstOcc, i * n, n);
                eulerLen[i] = tour.tourLen;
                leafCount[i] = tour.leafCount;
                bar.update(done.incrementAndGet());
            } catch (RuntimeException e) {
                throw treeBuildFailure(treeIndexBase + i, trees.get(i), e);
            }
        });
        bar.done();

        Config cfg = Config.getInstance();
        if (sm == null) sm = new SimilarityMatrix(n);
        GPUSimilarityMatrix.computeSimilarityGPUWide(
            depths, nodeSizes, true,
            null, null, null, null, null,
            micro, macro, firstOcc, eulerLen, leafCount,
            k, n, eMax, microLog, blockSize, blockMax, macroLog,
            cfg.getGpuDistTileSizeB(), effectiveTreeCapMiB(sm, cfg),
            cfg.getGpuDpProgressInterval(), cfg.getGpuDpProgressMaxSteps(),
            sm.numSum, sm.denSum,
            sm.isPacked() ? sm.packedNumeratorSegments() : null,
            sm.isPacked() ? sm.packedDenominatorSegments() : null,
            sm.packedSegmentShift());
        return sm;
    }

    public static SimilarityMatrix buildGPU(List<Tree> trees, int n) {
        int k = trees.size();
        preflightPackedHeap(n);
        if (SimilarityMatrix.requiresPacked(n)
                || Boolean.getBoolean("stelarpro.similarity.forcePacked")) {
            long triangle = SimilarityMatrix.triangleCellCount(n);
            Logging.info("Similarity output: exact segmented upper triangle "
                + "(%d cells, %.2f GiB per accumulator; dense n×n arrays are disabled)",
                triangle, triangle * Double.BYTES / (double)(1L << 30));
        }
        int eMaxRaw = 0;
        int logMax = 1;
        for (int i = 0; i < k; i++) {
            int len;
            try {
                len = EulerTourBuilder.tourLength(trees.get(i));
            } catch (RuntimeException e) {
                throw treeBuildFailure(i, trees.get(i), e);
            }
            eMaxRaw = Math.max(eMaxRaw, len);
            logMax = Math.max(logMax, sparseLog(len));
        }

        int ePadded = nextPowerOfTwo(eMaxRaw);
        boolean compactTourOverflow = eMaxRaw > COMPACT_MAX_TOUR;
        boolean compactLayoutOverflow = !fitsJavaArray((long)k * ePadded)
            || !fitsJavaArray((long)k * logMax * ePadded)
            || !fitsJavaArray((long)k * n);
        boolean forceWide = Boolean.getBoolean("stelarpro.similarity.forceWide");

        if (forceWide || compactTourOverflow || compactLayoutOverflow) {
            String reason = forceWide
                ? "forced by -Dstelarpro.similarity.forceWide=true"
                : compactTourOverflow
                    ? "maximum Euler tour " + eMaxRaw + " exceeds compact limit "
                        + COMPACT_MAX_TOUR
                    : "compact flat layout exceeds the Java single-array limit";
            Logging.info("GPU similarity RMQ: wide blocked mode (%s)", reason);
            return buildGPUWide(trees, n, eMaxRaw);
        }

        Logging.debug("GPU similarity RMQ: compact unsigned-16 mode (maximum tour %d)", eMaxRaw);
        return buildGPUCompact(trees, n);
    }

    private static SimilarityMatrix buildGPUCompact(List<Tree> trees, int n) {
        int k = trees.size();
        Logging.info("Building FullTourData for %d trees (parallel CPU)", k);

        // ── Step 1: Build per-tree full tour data in parallel ─────────────────
        EulerTourBuilder.FullTourData[] tours = new EulerTourBuilder.FullTourData[k];
        ProgressBar eulerBar    = new ProgressBar("FullTour + RMQ build", k);
        AtomicInteger eulerDone = new AtomicInteger(0);

        Threading.processRangeParallel(k, i -> {
            try {
                tours[i] = EulerTourBuilder.buildFull(trees.get(i), n);
                eulerBar.update(eulerDone.incrementAndGet());
            } catch (RuntimeException e) {
                throw treeBuildFailure(i, trees.get(i), e);
            }
        });
        eulerBar.done();

        // ── Step 2: Compute flat-array layout constants ───────────────────────
        int eMaxRaw = 0, logMaxRaw = 0;
        for (EulerTourBuilder.FullTourData td : tours) {
            if (td.tourLen > eMaxRaw)   eMaxRaw   = td.tourLen;
            if (td.log     > logMaxRaw) logMaxRaw = td.log;
        }
        int ePadded = 1;
        while (ePadded < eMaxRaw) ePadded <<= 1;
        final int E_max   = ePadded;
        final int LOG_max = logMaxRaw;

        // Per-position arrays:
        //   eulerDepths (short, 2B), eulerLeftChildS/RightChildS (short, 2B each),
        //   eulerF, eulerLeftChildF, eulerRightChildF (double, 8B each)
        // → 2 + 2 + 2 + 8 + 8 + 8 = 30 bytes/pos
        double euler_mb  = (double)k * E_max * 30 / 1e6;
        // Sparse table: one unsigned-16 left-biased argmin Euler position.
        // The kernel fetches depth and child payloads from the base Euler arrays.
        double sparse_mb = (double)k * LOG_max * E_max * Character.BYTES / 1e6;
        double leaf_mb   = (double)k * n * 4 / 1e6;
        Logging.info("  E_max=%d  LOG=%d  euler %.1f MB  sparse %.1f MB  leaf %.1f MB",
            E_max, LOG_max, euler_mb, sparse_mb, leaf_mb);

        // ── Step 3: Flatten into contiguous Java arrays ───────────────────────
        long edSize = (long)k * E_max;
        long spSize = (long)k * LOG_max * E_max;
        long ldSize = (long)k * n;

        short[]  eulerDepths       = new short [(int)edSize];
        double[] eulerF            = new double[(int)edSize];
        short[]  eulerLeftChildS   = new short [(int)edSize];
        double[] eulerLeftChildF   = new double[(int)edSize];
        short[]  eulerRightChildS  = new short [(int)edSize];
        double[] eulerRightChildF  = new double[(int)edSize];

        char[] sparseArgmin = new char[(int)spSize];

        int[]    firstOcc          = new int   [(int)ldSize];
        int[]    eulerLen          = new int   [k];
        int[]    leafCount         = new int   [k];

        Arrays.fill(firstOcc, -1);

        ProgressBar flatBar    = new ProgressBar("Flattening similarity tour data", k);
        AtomicInteger flatDone = new AtomicInteger(0);

        Threading.processRangeParallel(k, i -> {
            EulerTourBuilder.FullTourData td = tours[i];
            int len   = td.tourLen;
            int treeN = td.firstOcc.length;

            long edOff = (long)i * E_max;
            long spOff = (long)i * LOG_max * E_max;
            long ldOff = (long)i * n;

            for (int p = 0; p < len; p++) {
                int dst = (int)(edOff + p);
                eulerDepths      [dst] = td.depths[p];
                eulerF           [dst] = td.eulerF[p];
                eulerLeftChildS  [dst] = td.eulerLeftChildS[p];
                eulerLeftChildF  [dst] = td.eulerLeftChildF[p];
                eulerRightChildS [dst] = td.eulerRightChildS[p];
                eulerRightChildF [dst] = td.eulerRightChildF[p];
            }

            for (int lvl = 0; lvl < td.log; lvl++) {
                int rowLen = Math.max(0, len - (1 << lvl) + 1);
                long dst   = spOff + (long)lvl * E_max;
                for (int p = 0; p < rowLen; p++) {
                    int idx = (int)(dst + p);
                    sparseArgmin[idx] = td.sparseArgmin[lvl][p];
                }
            }

            for (int a = 0; a < treeN; a++) {
                int fo = td.firstOcc[a];
                if (fo >= 0) firstOcc[(int)(ldOff + a)] = fo;
            }

            eulerLen [i] = len;
            leafCount[i] = td.leafCount;
            flatBar.update(flatDone.incrementAndGet());
        });
        flatBar.done();

        // ── Step 4: Call GPU kernel ───────────────────────────────────────────
        Config cfg = Config.getInstance();
        int    tileSizeB        = cfg.getGpuDistTileSizeB();
        double progressInterval = cfg.getGpuDpProgressInterval();
        int    progressMaxSteps = cfg.getGpuDpProgressMaxSteps();

        SimilarityMatrix sm = new SimilarityMatrix(n);
        int treeCapMiB = effectiveTreeCapMiB(sm, cfg);
        GPUSimilarityMatrix.computeSimilarityGPU(
            eulerDepths,
            null, false,
            eulerF,
            eulerLeftChildS,  eulerLeftChildF,
            eulerRightChildS, eulerRightChildF,
            sparseArgmin,
            firstOcc, eulerLen, leafCount,
            k, n, E_max, LOG_max,
            tileSizeB, treeCapMiB,
            progressInterval, progressMaxSteps,
            sm.numSum, sm.denSum,
            sm.isPacked() ? sm.packedNumeratorSegments() : null,
            sm.isPacked() ? sm.packedDenominatorSegments() : null,
            sm.packedSegmentShift()
        );

        sm.normalize();
        return sm;
    }

    private static SimilarityMatrix buildGPUWide(List<Tree> trees, int n, int eMaxRaw) {
        if (!wideLayoutFits(trees.size(), n, eMaxRaw, MAX_JAVA_ARRAY_LENGTH)) {
            return buildGPUWideStreamed(trees, n, eMaxRaw);
        }
        SimilarityMatrix sm = accumulateGPUWide(trees, n, eMaxRaw, null, 0);
        sm.normalize();
        return sm;
    }

    /** Stream the established wide representation only when its one-shot arrays cannot fit. */
    private static SimilarityMatrix buildGPUWideStreamed(List<Tree> trees, int n,
                                                          int eMaxRaw) {
        int batchTrees = wideBatchTreeCount(trees.size(), n, eMaxRaw,
            MAX_JAVA_ARRAY_LENGTH, STREAMED_WIDE_HOST_BATCH_BYTES);
        int batches = (int)((trees.size() + (long)batchTrees - 1L) / batchTrees);
        Logging.info("Wide RMQ aggregate exceeds a Java array limit; streaming %d trees in "
                + "%d host batches (up to %d trees per batch; original order preserved)",
            trees.size(), batches, batchTrees);

        SimilarityMatrix sm = new SimilarityMatrix(n);
        sm.markStreamedHostBatches();
        for (int from = 0, batch = 1; from < trees.size(); batch++) {
            int to = (int)Math.min(trees.size(), from + (long)batchTrees);
            Logging.info("Wide similarity host batch %d/%d: trees %d..%d",
                batch, batches, from, to - 1);
            accumulateGPUWide(trees.subList(from, to), n, eMaxRaw, sm, from);
            from = to;
        }
        sm.normalize();
        return sm;
    }

    /** Flatten and accumulate one wide host batch; create an output matrix when {@code sm} is null. */
    private static SimilarityMatrix accumulateGPUWide(List<Tree> trees, int n, int eMaxRaw,
                                                       SimilarityMatrix sm,
                                                       int treeIndexBase) {
        int k = trees.size();
        int E_max = Math.max(1, eMaxRaw); // deliberately not power-of-two padded
        int blockSize = EulerTourBuilder.WIDE_BLOCK_SIZE;
        int microLog = EulerTourBuilder.WIDE_MICRO_LOG;
        int blockMax = (E_max + blockSize - 1) / blockSize;
        int macroLog = sparseLog(blockMax);

        int edSize = checkedLength((long)k * E_max, "wide Euler arrays");
        int microSize = checkedLength((long)k * microLog * E_max,
            "wide micro-RMQ array");
        int macroSize = checkedLength((long)k * macroLog * blockMax,
            "wide macro-RMQ array");
        int leafSize = checkedLength((long)k * n, "wide first-occurrence array");

        double baseGiB = (double)edSize * (Integer.BYTES * 3L + Double.BYTES * 3L)
            / (1L << 30);
        double rmqGiB = ((double)microSize + (double)macroSize * Integer.BYTES)
            / (1L << 30);
        Logging.info("Building wide FullTourData for %d trees (parallel CPU)", k);
        Logging.info("  E_max=%d  blocks=%d  microLOG=%d  macroLOG=%d  "
                + "Euler/payload %.2f GiB  RMQ %.2f GiB",
            E_max, blockMax, microLog, macroLog, baseGiB, rmqGiB);

        int[] eulerDepths = new int[edSize];
        double[] eulerF = new double[edSize];
        int[] eulerLeftChildS = new int[edSize];
        double[] eulerLeftChildF = new double[edSize];
        int[] eulerRightChildS = new int[edSize];
        double[] eulerRightChildF = new double[edSize];
        byte[] microArgmin = new byte[microSize];
        int[] macroArgmin = new int[macroSize];
        int[] firstOcc = new int[leafSize];
        int[] eulerLen = new int[k];
        int[] leafCount = new int[k];
        Arrays.fill(firstOcc, -1);

        ProgressBar bar = new ProgressBar("Wide FullTour + blocked RMQ build", k);
        AtomicInteger done = new AtomicInteger(0);
        Threading.processRangeParallel(k, i -> {
            Tree tree = trees.get(i);
            try {
                EulerTourBuilder.WideTourData td = EulerTourBuilder.buildWide(tree, n);
                int edOff = i * E_max;
                System.arraycopy(td.depths, 0, eulerDepths, edOff, td.tourLen);
                System.arraycopy(td.eulerF, 0, eulerF, edOff, td.tourLen);
                System.arraycopy(td.eulerLeftChildS, 0, eulerLeftChildS, edOff, td.tourLen);
                System.arraycopy(td.eulerLeftChildF, 0, eulerLeftChildF, edOff, td.tourLen);
                System.arraycopy(td.eulerRightChildS, 0, eulerRightChildS, edOff, td.tourLen);
                System.arraycopy(td.eulerRightChildF, 0, eulerRightChildF, edOff, td.tourLen);

                int microTreeOff = i * microLog * E_max;
                for (int lvl = 0; lvl < microLog; lvl++) {
                    System.arraycopy(td.microArgmin[lvl], 0, microArgmin,
                        microTreeOff + lvl * E_max, td.tourLen);
                }
                int macroTreeOff = i * macroLog * blockMax;
                for (int lvl = 0; lvl < td.macroLog; lvl++) {
                    int rowLen = Math.max(0, td.blockCount - (1 << lvl) + 1);
                    System.arraycopy(td.macroArgmin[lvl], 0, macroArgmin,
                        macroTreeOff + lvl * blockMax, rowLen);
                }
                System.arraycopy(td.firstOcc, 0, firstOcc, i * n, n);
                eulerLen[i] = td.tourLen;
                leafCount[i] = td.leafCount;
                bar.update(done.incrementAndGet());
            } catch (RuntimeException e) {
                throw treeBuildFailure(treeIndexBase + i, tree, e);
            }
        });
        bar.done();

        Config cfg = Config.getInstance();
        if (sm == null) sm = new SimilarityMatrix(n);
        int treeCapMiB = effectiveTreeCapMiB(sm, cfg);
        GPUSimilarityMatrix.computeSimilarityGPUWide(
            eulerDepths, null, false, eulerF,
            eulerLeftChildS, eulerLeftChildF,
            eulerRightChildS, eulerRightChildF,
            microArgmin, macroArgmin,
            firstOcc, eulerLen, leafCount,
            k, n, E_max, microLog, blockSize, blockMax, macroLog,
            cfg.getGpuDistTileSizeB(), treeCapMiB,
            cfg.getGpuDpProgressInterval(), cfg.getGpuDpProgressMaxSteps(),
            sm.numSum, sm.denSum,
            sm.isPacked() ? sm.packedNumeratorSegments() : null,
            sm.isPacked() ? sm.packedDenominatorSegments() : null,
            sm.packedSegmentShift());
        return sm;
    }

    static int effectiveTreeCapMiB(SimilarityMatrix sm, Config cfg) {
        int configured = cfg.getGpuSimilarityVramCapMiB();
        if (!sm.isPacked() || cfg.isGpuSimilarityVramCapExplicit()
                || configured >= LARGE_PACKED_GPU_TREE_CAP_MIB) return configured;
        Logging.info("Large-N similarity path: raising the tree-data batching ceiling "
            + "from %d MiB to %d MiB (still clamped to currently free VRAM)",
            configured, LARGE_PACKED_GPU_TREE_CAP_MIB);
        return LARGE_PACKED_GPU_TREE_CAP_MIB;
    }

    private static void preflightPackedHeap(int n) {
        if (!SimilarityMatrix.requiresPacked(n)
                && !Boolean.getBoolean("stelarpro.similarity.forcePacked")) return;
        long triangle = SimilarityMatrix.triangleCellCount(n);
        if (triangle > Long.MAX_VALUE / (2L * Double.BYTES)) {
            throw new IllegalArgumentException("Similarity matrix size overflows byte accounting");
        }
        long accumulatorBytes = triangle * 2L * Double.BYTES;
        Runtime rt = Runtime.getRuntime();
        long used = rt.totalMemory() - rt.freeMemory();
        long available = Math.max(0L, rt.maxMemory() - used);
        long safety = 256L << 20;
        if (available < accumulatorBytes + safety) {
            throw new IllegalArgumentException("Exact packed similarity accumulators for " + n
                + " taxa require at least "
                + String.format(java.util.Locale.ROOT, "%.2f", accumulatorBytes / (double)(1L << 30))
                + " GiB of currently available Java heap before similarity preprocessing; available "
                + String.format(java.util.Locale.ROOT, "%.2f", available / (double)(1L << 30))
                + " GiB. Increase --xmx or reduce the dataset size.");
        }
    }

    private static IllegalStateException treeBuildFailure(int index, Tree tree,
                                                           RuntimeException cause) {
        return new IllegalStateException("Similarity preprocessing failed for tree " + index
            + " (leaves=" + tree.leafCount + ", complete=" + tree.isComplete + ")", cause);
    }

    private static int sparseLog(int length) {
        if (length <= 1) return 1;
        return 32 - Integer.numberOfLeadingZeros(length);
    }

    private static int maximumTourLength(List<Tree> trees) {
        int maximum = 1;
        for (int i = 0; i < trees.size(); i++) {
            try {
                maximum = Math.max(maximum, EulerTourBuilder.tourLength(trees.get(i)));
            } catch (RuntimeException e) {
                throw treeBuildFailure(i, trees.get(i), e);
            }
        }
        return maximum;
    }

    private static int nextPowerOfTwo(int value) {
        if (value <= 1) return 1;
        if (value > (1 << 30)) return Integer.MAX_VALUE;
        return 1 << (32 - Integer.numberOfLeadingZeros(value - 1));
    }

    private static boolean fitsJavaArray(long length) {
        return length >= 0 && length <= MAX_JAVA_ARRAY_LENGTH;
    }

    /** Arithmetic-only check for the established flattened wide representation. */
    static boolean wideLayoutFits(int treeCount, int n, int eMaxRaw,
                                  long maxArrayLength) {
        if (treeCount < 0 || n <= 0 || eMaxRaw < 0 || maxArrayLength <= 0) {
            throw new IllegalArgumentException("invalid wide similarity layout dimensions");
        }
        int eMax = Math.max(1, eMaxRaw);
        int blockMax = (int)(((long)eMax + EulerTourBuilder.WIDE_BLOCK_SIZE - 1L)
            / EulerTourBuilder.WIDE_BLOCK_SIZE);
        int macroLog = sparseLog(blockMax);
        return (long)treeCount * eMax <= maxArrayLength
            && (long)treeCount * EulerTourBuilder.WIDE_MICRO_LOG * eMax <= maxArrayLength
            && (long)treeCount * macroLog * blockMax <= maxArrayLength
            && (long)treeCount * n <= maxArrayLength;
    }

    /**
     * Choose a safe streamed-wide batch size from Java's per-array element
     * bound and a total flattened-host-memory bound. Explicit limits keep the
     * boundary arithmetic testable without multi-gigabyte allocations.
     */
    static int wideBatchTreeCount(int treeCount, int n, int eMaxRaw,
                                  long maxArrayLength, long maxBatchBytes) {
        if (treeCount < 0 || n <= 0 || eMaxRaw < 0
                || maxArrayLength <= 0 || maxBatchBytes <= 0) {
            throw new IllegalArgumentException("invalid wide similarity batch dimensions");
        }
        if (treeCount == 0) return 1;

        int eMax = Math.max(1, eMaxRaw);
        int microLog = EulerTourBuilder.WIDE_MICRO_LOG;
        int blockMax = (int)(((long)eMax + EulerTourBuilder.WIDE_BLOCK_SIZE - 1L)
            / EulerTourBuilder.WIDE_BLOCK_SIZE);
        int macroLog = sparseLog(blockMax);
        long microCellsPerTree = Math.multiplyExact((long)microLog, eMax);
        long macroCellsPerTree = Math.multiplyExact((long)macroLog, blockMax);
        long flattenedBytesPerTree = Math.addExact(
            Math.addExact(Math.multiplyExact((long)eMax, 36L), microCellsPerTree),
            Math.addExact(Math.multiplyExact(macroCellsPerTree, Integer.BYTES),
                Math.addExact(Math.multiplyExact((long)n, Integer.BYTES),
                    2L * Integer.BYTES)));

        long byEulerArray = maxArrayLength / eMax;
        long byMicroArray = maxArrayLength / microCellsPerTree;
        long byMacroArray = maxArrayLength / macroCellsPerTree;
        long byFirstOccurrence = maxArrayLength / n;
        long byHostBytes = maxBatchBytes / flattenedBytesPerTree;
        long safe = Math.min(Math.min(byEulerArray, byMicroArray),
            Math.min(Math.min(byMacroArray, byFirstOccurrence), byHostBytes));
        if (safe < 1) {
            throw new IllegalArgumentException("one wide similarity tree cannot fit in the "
                + "configured streamed host-array bounds");
        }
        return (int)Math.min(treeCount, Math.min(safe, Integer.MAX_VALUE));
    }

    private static int checkedLength(long length, String label) {
        if (!fitsJavaArray(length)) {
            throw new IllegalArgumentException(label + " requires " + length
                + " elements; Java arrays support at most " + MAX_JAVA_ARRAY_LENGTH
                + ". Reduce the number of trees per run or use a future streamed layout.");
        }
        return (int)length;
    }

    // ── CPU accumulation (per-tree, scatter mirroring ASTRAL-MP) ─────────────

    /**
     * Scatter same-side quartet counts to all (l, r) pairs at every internal
     * node, plus accumulate the per-pair denominator C2(kt − 2).
     */
    private static void accumulateCPU(Tree tree, SimilarityMatrix sm) {
        int kt = tree.leafCount;
        if (kt < 4) {
            // C2(kt-2) = 0 → tree contributes nothing to numerator or denominator
            // (matches ASTRAL-MP, which skips by virtue of sim and den both being 0).
            return;
        }

        long denPerPair = EulerTourBuilder.c2(kt - 2);
        int n = sm.n;
        int[] postArr = tree.postorderArray;

        // ── Numerator: scatter across all internal nodes ─────────────────────
        if (sm.isPacked()) scatterAtNodePacked(tree.root, tree, kt, sm);
        else               scatterAtNodeDense(tree.root, tree, kt, sm);

        // ── Denominator: every pair (a, b) co-occurring in T gets += C2(kt-2) ─
        // ASTRAL-MP accumulates 2·C2(kt-2) in dn[l][r] (doubled by the
        // mirror-write pattern), then normalizes by dn/2. Equivalent to a
        // single sum here; the constant factor 2 cancels.
        if (sm.isPacked()) {
            for (int i = 0; i < kt; i++) {
                int a = postArr[i];
                for (int j = i + 1; j < kt; j++) {
                    sm.addPackedDenominator(a, postArr[j], denPerPair);
                }
            }
        } else {
            for (int i = 0; i < kt; i++) {
                int a = postArr[i];
                long rowOff = (long) a * n;
                for (int j = 0; j < kt; j++) {
                    if (i == j) continue;
                    int b = postArr[j];
                    sm.denSum[(int)(rowOff + b)] += denPerPair;
                }
            }
        }
    }

    /**
     * Post-order scatter at one internal node.
     *
     * Components at u:
     *   left   = subtree of u.left   (rangeStart_L .. rangeEnd_L)
     *   right  = subtree of u.right  (rangeStart_R .. rangeEnd_R)
     *   others = leaves of T outside subtree(u)
     *            (i.e. positions [0, u.rangeStart) ∪ [u.rangeEnd, kt) in postArr)
     *
     * For each ordered component-pair, scatter
     *   sim = totalPairs − C2(|X|) − C2(|Y|)
     * to every (l ∈ X, r ∈ Y) leaf pair, symmetrically.
     */
    private static void scatterAtNodeDense(TreeNode node, Tree tree, int kt,
                                           SimilarityMatrix sm) {
        if (node.isLeaf()) return;
        scatterAtNodeDense(node.left,  tree, kt, sm);
        scatterAtNodeDense(node.right, tree, kt, sm);

        int subL = node.left.rangeEnd  - node.left.rangeStart;
        int subR = node.right.rangeEnd - node.right.rangeStart;
        int subU = node.rangeEnd       - node.rangeStart;
        int subO = kt - subU;                                 // "others" size

        long cL = EulerTourBuilder.c2(subL);
        long cR = EulerTourBuilder.c2(subR);
        long cO = EulerTourBuilder.c2(subO);
        long totalPairs = cL + cR + cO;                       // only positive comps survive

        int n = sm.n;
        int[] postArr = tree.postorderArray;

        // ── (left × right): always present ───────────────────────────────────
        long simLR = totalPairs - cL - cR;
        if (simLR != 0) {
            scatterRangeRange(postArr,
                node.left.rangeStart, node.left.rangeEnd,
                node.right.rangeStart, node.right.rangeEnd,
                simLR, sm.numSum, n);
        }

        // ── (left × others) and (right × others): only when u is non-root ────
        if (subO > 0) {
            long simLO = totalPairs - cL - cO;
            if (simLO != 0) {
                scatterRangeOthers(postArr, kt,
                    node.left.rangeStart, node.left.rangeEnd,
                    node.rangeStart, node.rangeEnd, simLO, sm.numSum, n);
            }
            long simRO = totalPairs - cR - cO;
            if (simRO != 0) {
                scatterRangeOthers(postArr, kt,
                    node.right.rangeStart, node.right.rangeEnd,
                    node.rangeStart, node.rangeEnd, simRO, sm.numSum, n);
            }
        }
    }

    private static void scatterAtNodePacked(TreeNode node, Tree tree, int kt,
                                            SimilarityMatrix sm) {
        if (node.isLeaf()) return;
        scatterAtNodePacked(node.left, tree, kt, sm);
        scatterAtNodePacked(node.right, tree, kt, sm);

        int subL = node.left.rangeEnd - node.left.rangeStart;
        int subR = node.right.rangeEnd - node.right.rangeStart;
        int subU = node.rangeEnd - node.rangeStart;
        int subO = kt - subU;
        long cL = EulerTourBuilder.c2(subL);
        long cR = EulerTourBuilder.c2(subR);
        long cO = EulerTourBuilder.c2(subO);
        long totalPairs = cL + cR + cO;
        int[] postArr = tree.postorderArray;

        long simLR = totalPairs - cL - cR;
        if (simLR != 0) {
            scatterRangeRangePacked(postArr,
                node.left.rangeStart, node.left.rangeEnd,
                node.right.rangeStart, node.right.rangeEnd, simLR, sm);
        }
        if (subO > 0) {
            long simLO = totalPairs - cL - cO;
            if (simLO != 0) {
                scatterRangeOthersPacked(postArr, kt,
                    node.left.rangeStart, node.left.rangeEnd,
                    node.rangeStart, node.rangeEnd, simLO, sm);
            }
            long simRO = totalPairs - cR - cO;
            if (simRO != 0) {
                scatterRangeOthersPacked(postArr, kt,
                    node.right.rangeStart, node.right.rangeEnd,
                    node.rangeStart, node.rangeEnd, simRO, sm);
            }
        }
    }

    /** Scatter `sim` to every (a ∈ [aLo,aHi)) × (b ∈ [bLo,bHi)) leaf-pair, both directions. */
    private static void scatterRangeRange(int[] postArr,
                                           int aLo, int aHi,
                                           int bLo, int bHi,
                                           long sim, double[] numSum, int n) {
        double sd = (double) sim;
        for (int pi = aLo; pi < aHi; pi++) {
            int a = postArr[pi];
            long rowA = (long) a * n;
            for (int pj = bLo; pj < bHi; pj++) {
                int b = postArr[pj];
                numSum[(int)(rowA + b)] += sd;
                numSum[b * n + a]       += sd;
            }
        }
    }

    private static void scatterRangeRangePacked(int[] postArr,
                                                 int aLo, int aHi,
                                                 int bLo, int bHi,
                                                 long sim, SimilarityMatrix sm) {
        double sd = (double)sim;
        for (int pi = aLo; pi < aHi; pi++) {
            int a = postArr[pi];
            for (int pj = bLo; pj < bHi; pj++) {
                sm.addPackedNumerator(a, postArr[pj], sd);
            }
        }
    }

    /**
     * Scatter `sim` to every (a ∈ [aLo,aHi)) × (b ∈ others) leaf-pair, where
     * others = leaves in T outside [subLo, subHi).
     */
    private static void scatterRangeOthers(int[] postArr, int kt,
                                            int aLo, int aHi,
                                            int subLo, int subHi,
                                            long sim, double[] numSum, int n) {
        double sd = (double) sim;
        for (int pi = aLo; pi < aHi; pi++) {
            int a = postArr[pi];
            long rowA = (long) a * n;
            for (int pj = 0; pj < subLo; pj++) {
                int b = postArr[pj];
                numSum[(int)(rowA + b)] += sd;
                numSum[b * n + a]       += sd;
            }
            for (int pj = subHi; pj < kt; pj++) {
                int b = postArr[pj];
                numSum[(int)(rowA + b)] += sd;
                numSum[b * n + a]       += sd;
            }
        }
    }

    private static void scatterRangeOthersPacked(int[] postArr, int kt,
                                                  int aLo, int aHi,
                                                  int subLo, int subHi,
                                                  long sim, SimilarityMatrix sm) {
        double sd = (double)sim;
        for (int pi = aLo; pi < aHi; pi++) {
            int a = postArr[pi];
            for (int pj = 0; pj < subLo; pj++) {
                sm.addPackedNumerator(a, postArr[pj], sd);
            }
            for (int pj = subHi; pj < kt; pj++) {
                sm.addPackedNumerator(a, postArr[pj], sd);
            }
        }
    }
}
