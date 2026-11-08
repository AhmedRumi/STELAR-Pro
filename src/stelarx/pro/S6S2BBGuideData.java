package stelarx.pro;

import stelarx.Config;
import stelarx.Logging;
import stelarx.completion.SimilarityMatrix;
import stelarx.completion.SimilarityMatrixBuilder;
import stelarx.completion.SortedRows;
import stelarx.completion.SortedRowsBuilder;
import stelarx.completion.UPGMAClusterer;
import stelarx.gpu.GPUSimilarityMatrix;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Shared strict-lineage S4 guide pool for every S6 backbone. */
public final class S6S2BBGuideData {
    private final SimilarityMatrix matrix;
    private final List<Tree> completed;
    private final Tree upgma;
    private final TaxonRegistry global;

    private S6S2BBGuideData(SimilarityMatrix matrix, List<Tree> completed,
                            Tree upgma, TaxonRegistry global) {
        this.matrix = matrix;
        this.completed = completed;
        this.upgma = upgma;
        this.global = global;
    }

    public static S6S2BBGuideData prepare(List<Tree> genes, TaxonRegistry taxa,
                                          int upgmaIndex) {
        long start = System.nanoTime();
        List<Tree> disco = DiscoDecomposer.decomposeAll(
            genes, taxa.size(), 0, 4, true).trees();
        if (disco.isEmpty())
            throw new IllegalStateException("S6-S2BB needs at least one DISCO component");
        boolean gpu = Config.getInstance().getComputeMode() == Config.ComputeMode.GPU
            && GPUSimilarityMatrix.tryLoad();
        // S4 uses rooted-triplet sister frequencies for both completion and
        // UPGMA; reuse its bounded JNI/GPU host batches for this shared pool.
        SimilarityMatrix matrix = gpu
            ? SimilarityMatrixBuilder.buildTripletGPU(disco, taxa.size())
            : SimilarityMatrixBuilder.buildTripletCPU(disco, taxa.size());
        SortedRows rows = ScopedTreeCompleter.anyMissing(disco)
            ? SortedRowsBuilder.build(matrix) : null;
        List<Tree> completed = ScopedTreeCompleter.prepare(rows, disco, taxa.size())
            .completeBatch(disco);
        // UPGMA receives a private triangle; the packed global matrix remains
        // read-only for group-sized UPGMA partition trees in all four rounds.
        Tree upgma = UPGMAClusterer.build(matrix, upgmaIndex, false);
        markCandidate(upgma.root);
        Logging.info("S6-S2BB shared S4 guide pool: %d strict-lineage DISCO trees, "
            + "%d taxa, %s rooted-triplet similarity, one UPGMA tree in %d ms",
            completed.size(), taxa.size(), gpu ? "GPU" : "CPU",
            (System.nanoTime() - start) / 1_000_000);
        return new S6S2BBGuideData(matrix, completed, upgma, taxa);
    }

    public SimilarityMatrix matrix() { return matrix; }
    public Tree upgmaTree() { return upgma; }

    /** Induce completed DISCO/UPGMA proposals on one random backbone subset. */
    public List<Tree> restrict(TaxonRegistry subset) {
        int[] map = new int[global.size()];
        Arrays.fill(map, -1);
        for (int i = 0; i < subset.size(); i++)
            map[global.getId(subset.getName(i))] = i;
        List<Tree> guides = new ArrayList<>(completed.size() + 1);
        for (Tree source : completed) {
            Tree induced = restrictOne(source, map, guides.size(), subset.size());
            if (induced != null && induced.leafCount >= 3) guides.add(induced);
        }
        Tree inducedUpgma = restrictOne(upgma, map, guides.size(), subset.size());
        if (inducedUpgma == null || !inducedUpgma.isComplete)
            throw new IllegalStateException("S6-S2BB subset UPGMA is incomplete");
        guides.add(inducedUpgma);
        Logging.info("S6-S2BB backbone: %d taxa, %d restricted DISCO + 1 UPGMA proposals",
            subset.size(), guides.size() - 1);
        return guides;
    }

    private static Tree restrictOne(Tree source, int[] map, int index, int n) {
        TreeNode root = prune(source.root, map);
        if (root == null) return null;
        int[] order = new int[n];
        int[] count = {0};
        assignRanges(root, order, count);
        int[] leaves = Arrays.copyOf(order, count[0]);
        int[] position = new int[n];
        Arrays.fill(position, -1);
        for (int i = 0; i < leaves.length; i++) {
            if (position[leaves[i]] != -1)
                throw new IllegalStateException("S6-S2BB subset duplicated a taxon");
            position[leaves[i]] = i;
        }
        return new Tree(index, root, leaves, position, leaves.length, n);
    }

    private static TreeNode prune(TreeNode node, int[] map) {
        if (node.isLeaf()) {
            int id = map[node.taxonId];
            if (id < 0) return null;
            TreeNode leaf = new TreeNode();
            leaf.taxonId = id;
            return leaf;
        }
        if (node.isPolytomous())
            throw new IllegalStateException("S6-S2BB requires binary completed DISCO trees");
        TreeNode left = prune(node.left, map);
        TreeNode right = prune(node.right, map);
        if (left == null) return right;
        if (right == null) return left;
        TreeNode parent = new TreeNode();
        parent.left = left;
        parent.right = right;
        left.parent = right.parent = parent;
        parent.isSpeciationNode = true;
        parent.isDuplicationNode = false;
        return parent;
    }

    private static void assignRanges(TreeNode node, int[] order, int[] cursor) {
        node.rangeStart = cursor[0];
        if (node.isLeaf()) order[cursor[0]++] = node.taxonId;
        else {
            assignRanges(node.left, order, cursor);
            assignRanges(node.right, order, cursor);
        }
        node.rangeEnd = cursor[0];
    }

    private static void markCandidate(TreeNode node) {
        if (node.isLeaf()) return;
        node.isSpeciationNode = true;
        node.isDuplicationNode = false;
        markCandidate(node.left);
        markCandidate(node.right);
    }
}
