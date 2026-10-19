package stelarx.pro;

import stelarx.Logging;
import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.completion.SimilarityMatrix;
import stelarx.completion.TreeCompleter;
import stelarx.dp.DPTable;
import stelarx.hash.TaxonHasher;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/** Streams S2's completed DISCO proposals into X and the rooted DP space. */
public final class S2CandidateStreamer {
    private S2CandidateStreamer() {}

    public record Result(List<Tree> exemplarTrees, int newClusters,
                         int completedTrees, int retainedExemplars, long transitions) {}

    /**
     * Completion remains parallel within each bounded batch. Extraction is the
     * same postorder single-copy hash/transition walk previously performed in
     * Phases 3 and 5; only leaf indexes of trees that introduce new clusters
     * survive for later CPU/CUDA membership lookups.
     */
    public static Result add(List<Tree> discoTrees, SimilarityMatrix similarity,
                             Tree upgma, List<Tree> originalTrees,
                             ClusterTable clusters, DPTable dp,
                             TaxonHasher hasher, TaxonRegistry registry,
                             int batchSize, String completedTreesFile) throws IOException {
        if (batchSize < 1) throw new IllegalArgumentException("S2 batch size must be positive");
        int n = registry.size();
        if (discoTrees.isEmpty() && upgma == null && completedTreesFile == null) {
            return new Result(originalTrees, 0, 0, 0, 0);
        }
        List<Tree> exemplars = new ArrayList<>(originalTrees.size() + discoTrees.size() + 1);
        exemplars.addAll(originalTrees);
        ClusterHash[] leaves = singletonHashes(hasher);
        TreeCompleter.PreparedCompletion completer = discoTrees.isEmpty()
            ? null : TreeCompleter.prepare(similarity, n);
        int initialClusters = clusters.size();
        int newClusters = 0;
        long transitions = 0;
        int retained = 0;

        try (PrintStream dump = completedTreesFile == null ? null
                : new PrintStream(new FileOutputStream(completedTreesFile))) {
            for (int start = 0; start < discoTrees.size(); start += batchSize) {
                int end = Math.min(start + batchSize, discoTrees.size());
                List<Tree> completed = completer.completeBatch(discoTrees.subList(start, end));
                for (Tree tree : completed) {
                    if (!tree.isComplete || tree.leafCount != n) {
                        throw new IllegalStateException("S2 completion did not yield a complete "
                            + "single-copy tree " + tree.treeIndex);
                    }
                    if (dump != null) dump.println(tree.toNewick(registry));
                    int index = exemplars.size();
                    Counts counts = new Counts();
                    ClusterHash root = emit(tree.root, tree, index, leaves,
                        hasher.numSeeds(), clusters, dp, counts);
                    if (!root.equals(clusters.getAllTaxaHash())) {
                        throw new IllegalStateException("S2 completed tree has an invalid root hash");
                    }
                    newClusters += counts.newClusters;
                    transitions += counts.transitions;
                    if (counts.newClusters > 0) {
                        exemplars.add(Tree.candidateExemplar(index, tree, n));
                        retained++;
                    }
                }
                // Drop the uncompleted DISCO roots too; the mutable list keeps
                // only not-yet-processed components alive.
                for (int i = start; i < end; i++) discoTrees.set(i, null);
                // No completed node graph or per-node hash index survives this batch.
                completed = null;
                if (end == discoTrees.size() || (end / batchSize) % 100 == 0) {
                    // UPGMA is emitted after the DISCO batches, so this is not
                    // the final S2 candidate total.
                    Logging.info("S2 DISCO progress: %d/%d trees, +%d unique clusters "
                            + "(UPGMA pending)",
                        end, discoTrees.size(), newClusters);
                }
            }
        }

        int discoNewClusters = newClusters;
        int upgmaNewClusters = 0;
        if (upgma != null) {
            Counts counts = new Counts();
            int index = exemplars.size();
            emit(upgma.root, upgma, index, leaves, hasher.numSeeds(), clusters, dp, counts);
            newClusters += counts.newClusters;
            upgmaNewClusters = counts.newClusters;
            transitions += counts.transitions;
            if (counts.newClusters > 0) {
                exemplars.add(Tree.candidateExemplar(index, upgma, n));
                retained++;
            }
        }
        if (completedTreesFile != null) {
            Logging.info("Completed gene trees written to %s (%d trees)",
                completedTreesFile, discoTrees.size());
        }
        if (clusters.size() != initialClusters + newClusters) {
            throw new IllegalStateException("S2 streamed candidate count disagrees with X");
        }
        Logging.info("S2 streamed candidates: %d DISCO + %d UPGMA proposal(s), "
                + "+%d unique clusters (%d DISCO, %d UPGMA; %d total in X), "
                + "%d local transitions, %d compact exemplars retained",
            discoTrees.size(), upgma == null ? 0 : 1,
            newClusters, discoNewClusters, upgmaNewClusters, clusters.size(),
            transitions, retained);
        return new Result(exemplars, newClusters, discoTrees.size(), retained, transitions);
    }

    private static ClusterHash[] singletonHashes(TaxonHasher hasher) {
        int n = hasher.numTaxa(), m = hasher.numSeeds();
        ClusterHash[] hashes = new ClusterHash[n];
        for (int taxon = 0; taxon < n; taxon++) {
            long[] sums = new long[m], xors = new long[m];
            for (int seed = 0; seed < m; seed++) {
                sums[seed] = xors[seed] = hasher.get(seed, taxon);
            }
            hashes[taxon] = new ClusterHash(sums, xors, 1, m);
        }
        return hashes;
    }

    private static final class Counts {
        int newClusters;
        long transitions;
    }

    private static ClusterHash emit(TreeNode node, Tree tree, int exemplarIndex,
                                    ClusterHash[] leaves, int seeds,
                                    ClusterTable clusters, DPTable dp, Counts counts) {
        ClusterHash hash;
        if (node.isLeaf()) {
            hash = leaves[node.taxonId];
        } else {
            if (node.isPolytomous()) {
                throw new IllegalStateException("Completed S2 proposal is not binary");
            }
            ClusterHash left = emit(node.left, tree, exemplarIndex,
                leaves, seeds, clusters, dp, counts);
            ClusterHash right = emit(node.right, tree, exemplarIndex,
                leaves, seeds, clusters, dp, counts);
            long[] sums = new long[seeds], xors = new long[seeds];
            for (int seed = 0; seed < seeds; seed++) {
                sums[seed] = left.sums[seed] + right.sums[seed];
                xors[seed] = left.xors[seed] ^ right.xors[seed];
            }
            hash = new ClusterHash(sums, xors, left.size + right.size, seeds);
            dp.addStreamedTransition(hash, left, right);
            counts.transitions++;
        }
        if (!node.isRoot() && clusters.addStreamedCandidate(hash,
                exemplarIndex, node.rangeStart, node.rangeEnd, tree.leafCount)) {
            counts.newClusters++;
        }
        return hash;
    }
}
