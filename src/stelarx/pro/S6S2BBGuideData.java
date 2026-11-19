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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** S6 subset-local guides with duplication-vertex completion scopes. */
public final class S6S2BBGuideData {
    private final TaxonRegistry global;
    private final List<Tree> uncompleted;

    private S6S2BBGuideData(TaxonRegistry global,
                            List<Tree> uncompleted) {
        this.global = global;
        this.uncompleted = uncompleted;
    }

    public static S6S2BBGuideData prepare(List<Tree> genes, TaxonRegistry taxa) {
        List<Tree> disco = DiscoDecomposer.decomposeAll(
            genes, taxa.size(), 0, 4, DiscoDecomposer.CompletionScope.DUPLICATION).trees();
        if (disco.isEmpty())
            throw new IllegalStateException("S6-S2BB needs at least one DISCO component");
        Logging.info("S6-local: retaining %d uncompleted DISCO components; "
            + "duplication-vertex scope; no global similarity, completion, or UPGMA",
            disco.size());
        return new S6S2BBGuideData(taxa, disco);
    }

    /** Restrict first, then complete and build UPGMA only inside this subset. */
    public List<Tree> restrict(TaxonRegistry subset) {
        return restrictCompleteLocally(subset);
    }

    /**
     * Restrict uncompleted components before computing the subset's triplet
     * similarity. Missing species are intersected with the pre-cut duplication
     * scope; taxa above that duplication are excluded. Scoping affects candidate
     * completion only; original tagged genes supply all triplet weights.
     */
    private List<Tree> restrictCompleteLocally(TaxonRegistry subset) {
        int n = subset.size();
        int[] map = new int[global.size()];
        Arrays.fill(map, -1);
        for (int i = 0; i < n; i++)
            map[global.getId(subset.getName(i))] = i;
        List<Tree> components = new ArrayList<>();
        List<Tree> observed = new ArrayList<>();
        for (Tree source : uncompleted) {
            TreeNode root = prune(source.root, map);
            if (root == null) continue; // No present anchor in this subset.
            long[] missing = restrictMissing((DiscoDecomposer.ScopedTree) source,
                map, n);
            DiscoDecomposer.ScopedTree scoped =
                DiscoDecomposer.restrictedScopedTree(components.size(), root, n, missing);
            if (scoped.leafCount + scoped.missingCount() < 3) continue;
            components.add(scoped);
            if (scoped.leafCount >= 3) observed.add(scoped);
        }
        SimilarityMatrix local = tripletMatrix(observed, n, true);
        SortedRows rows = ScopedTreeCompleter.anyMissing(components)
            ? SortedRowsBuilder.build(local) : null;
        List<Tree> completedHere = ScopedTreeCompleter.prepare(rows, components, n)
            .completeBatch(components);
        List<Tree> guides = new ArrayList<>(completedHere.size() + 1);
        for (Tree tree : completedHere)
            if (tree.leafCount >= 3) guides.add(tree);
        Tree localUpgma = UPGMAClusterer.build(local, guides.size(), true);
        markCandidate(localUpgma.root);
        guides.add(localUpgma);
        Logging.info("S6-local backbone: %d taxa, %d restricted completed DISCO "
            + "+ 1 local UPGMA proposals", n, guides.size() - 1);
        return guides;
    }

    private static long[] restrictMissing(DiscoDecomposer.ScopedTree source,
                                          int[] map, int n) {
        long[] result = new long[(n + 63) >>> 6];
        long[] sourceWords = source.missingWords();
        for (int wi = 0; wi < sourceWords.length; wi++) {
            long word = sourceWords[wi];
            while (word != 0) {
                int original = (wi << 6) + Long.numberOfTrailingZeros(word);
                word &= word - 1;
                if (original < map.length && map[original] >= 0) {
                    int local = map[original];
                    result[local >>> 6] |= 1L << local;
                }
            }
        }
        return result;
    }

    private static SimilarityMatrix tripletMatrix(List<Tree> observed, int n,
                                                   boolean allowGpu) {
        if (observed.isEmpty()) {
            SimilarityMatrix empty = new SimilarityMatrix(n);
            empty.normalize();
            return empty;
        }
        // Tiny partition matrices are cheaper on CPU than one CUDA launch each.
        boolean gpu = allowGpu && n >= 64
            && Config.getInstance().getComputeMode() == Config.ComputeMode.GPU
            && GPUSimilarityMatrix.tryLoad();
        return gpu ? SimilarityMatrixBuilder.buildTripletGPU(observed, n)
            : SimilarityMatrixBuilder.buildTripletCPU(observed, n);
    }

    /**
     * After edge assignment, all partition groups are disjoint. Route every
     * observed DISCO leaf to its group in one postorder traversal per component
     * (small-to-large maps), then build one small matrix per nontrivial group.
     */
    public SimilarityMatrix[] partitionMatrices(List<List<Integer>> groups) {
        int count = groups.size();
        int[] groupOf = new int[global.size()];
        int[] localOf = new int[global.size()];
        Arrays.fill(groupOf, -1);
        for (int g = 0; g < count; g++) {
            List<Integer> taxa = groups.get(g);
            if (taxa.size() < 3) continue;
            for (int j = 0; j < taxa.size(); j++) {
                int taxon = taxa.get(j);
                if (groupOf[taxon] >= 0)
                    throw new IllegalStateException("Overlapping S6 partition groups");
                groupOf[taxon] = g;
                localOf[taxon] = j;
            }
        }
        List<List<Tree>> sources = new ArrayList<>(count);
        for (int g = 0; g < count; g++) sources.add(new ArrayList<>());
        for (Tree source : uncompleted) {
            Map<Integer, TreeNode> roots = distribute(source.root, groupOf, localOf);
            for (Map.Entry<Integer, TreeNode> entry : roots.entrySet()) {
                int g = entry.getKey(), size = groups.get(g).size();
                List<Tree> target = sources.get(g);
                Tree restricted = DiscoDecomposer.restrictedScopedTree(
                    target.size(), entry.getValue(), size, new long[0]);
                if (restricted.leafCount >= 3) target.add(restricted);
            }
        }
        SimilarityMatrix[] matrices = new SimilarityMatrix[count];
        for (int g = 0; g < count; g++) {
            if (groups.get(g).size() >= 3)
                matrices[g] = tripletMatrix(sources.get(g), groups.get(g).size(), false);
            sources.set(g, null); // release each group's restricted tree nodes
        }
        Logging.info("S6-local partition matrices: %d disjoint groups from one "
            + "DISCO restriction pass", count);
        return matrices;
    }

    private static Map<Integer, TreeNode> distribute(TreeNode node,
                                                      int[] groupOf, int[] localOf) {
        if (node.isLeaf()) {
            int taxon = node.taxonId;
            int group = groupOf[taxon];
            if (group < 0) return new HashMap<>(0);
            TreeNode leaf = new TreeNode();
            leaf.taxonId = localOf[taxon];
            Map<Integer, TreeNode> one = new HashMap<>(1);
            one.put(group, leaf);
            return one;
        }
        Map<Integer, TreeNode> left = distribute(node.left, groupOf, localOf);
        Map<Integer, TreeNode> right = distribute(node.right, groupOf, localOf);
        boolean leftIsLarge = left.size() >= right.size();
        Map<Integer, TreeNode> large = leftIsLarge ? left : right;
        Map<Integer, TreeNode> small = leftIsLarge ? right : left;
        for (Map.Entry<Integer, TreeNode> entry : small.entrySet()) {
            int group = entry.getKey();
            TreeNode old = large.get(group);
            if (old == null) {
                large.put(group, entry.getValue());
            } else {
                TreeNode parent = new TreeNode();
                parent.left = leftIsLarge ? old : entry.getValue();
                parent.right = leftIsLarge ? entry.getValue() : old;
                parent.left.parent = parent.right.parent = parent;
                parent.isSpeciationNode = true;
                large.put(group, parent);
            }
        }
        return large;
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
            throw new IllegalStateException("S6-S2BB requires binary DISCO components");
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

    private static void markCandidate(TreeNode node) {
        if (node.isLeaf()) return;
        node.isSpeciationNode = true;
        node.isDuplicationNode = false;
        markCandidate(node.left);
        markCandidate(node.right);
    }
}
