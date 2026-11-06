package stelarx.completion;

import stelarx.Logging;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.util.ProgressBar;
import stelarx.util.Threading;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Inserts missing taxa into incomplete rooted gene trees while preserving the
 * supplied root and every relationship among the taxa already present.
 *
 * Algorithm for inserting taxon x into tree T:
 *   1. Find anchor a = closest taxon currently present in T (via {@link SortedRows}).
 *   2. Replace leaf a by a new binary parent whose children are a and x.
 *   3. Update the present-taxon and leaf-node maps.
 *
 * This never reverses an edge. Deleting inserted leaves and suppressing their
 * degree-one parents therefore recovers exactly the supplied rooted tree.
 *
 * Missing taxa are processed in ascending ID order for each tree.
 * Different trees are completed in parallel using the thread pool.
 *
 * The only matrix-derived input is the nearest-neighbour order of each taxon
 * ({@link SortedRows}); similarity values themselves are never read here, so
 * the similarity triangle can be released (or consumed by UPGMA) once the rows
 * exist.
 */
public class TreeCompleter {

    /** Reuses the read-only nearest-neighbour index across bounded completion batches. */
    public static final class PreparedCompletion {
        private final SortedRows rows;
        private final int numTaxa;

        private PreparedCompletion(SortedRows rows, int numTaxa) {
            this.rows = rows;
            this.numTaxa = numTaxa;
        }

        /** Complete independent trees in parallel; only this batch is retained. */
        public List<Tree> completeBatch(List<Tree> trees) {
            Tree[] result = trees.toArray(new Tree[0]);
            Threading.processRangeParallel(result.length, i -> {
                if (!result[i].isComplete) {
                    result[i] = completeTreeRootPreserving(result[i], rows, numTaxa);
                }
            });
            return Arrays.asList(result);
        }
    }

    /** Prepare completion over previously built nearest-neighbour rows. */
    public static PreparedCompletion prepare(SortedRows rows, int numTaxa) {
        if (rows == null) throw new IllegalArgumentException("completion rows are required");
        if (rows.n != numTaxa) {
            throw new IllegalArgumentException("completion rows cover " + rows.n
                + " taxa, expected " + numTaxa);
        }
        return new PreparedCompletion(rows, numTaxa);
    }

    /** True when at least one tree lacks a taxon. */
    public static boolean anyIncomplete(List<Tree> trees) {
        for (Tree tree : trees) if (!tree.isComplete) return true;
        return false;
    }

    /**
     * Complete all incomplete trees with root-preserving nearest-neighbor insertion.
     *
     * @param trees      gene trees (may mix complete and incomplete)
     * @param rows       nearest-neighbour order per taxon (see {@link SortedRowsBuilder})
     * @param n          total taxon count
     * @return new list where every tree is complete; already-complete trees pass through
     */
    public static List<Tree> completeAll(List<Tree> trees, SortedRows rows, int n) {
        List<Integer> incomplete = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            if (!trees.get(i).isComplete) incomplete.add(i);
        }
        if (incomplete.isEmpty()) return trees;
        if (rows == null) throw new IllegalArgumentException("completion rows are required");

        Logging.info("Tree completion: %d/%d trees incomplete", incomplete.size(), trees.size());

        // Mutable result array; complete trees pass through unchanged.
        Tree[] result = trees.toArray(new Tree[0]);
        ProgressBar bar   = new ProgressBar("Completing incomplete gene trees", incomplete.size());
        AtomicInteger cnt = new AtomicInteger(0);
        Threading.processParallel(incomplete, idx -> {
            result[idx] = completeTreeRootPreserving(trees.get(idx), rows, n);
            bar.update(cnt.incrementAndGet());
        });
        bar.done();
        return Arrays.asList(result);
    }

    /**
     * Convenience entry point: builds the nearest-neighbour rows from the
     * matrix only when some tree is incomplete.  The matrix is only read.
     */
    public static List<Tree> completeAll(List<Tree> trees, SimilarityMatrix sim, int n) {
        if (!anyIncomplete(trees)) return trees;
        return completeAll(trees, SortedRowsBuilder.build(sim), n);
    }

    /**
     * Root-aware completion used by STELAR-Pro. A missing taxon is attached as the
     * sibling of its most similar present taxon. Replacing that leaf by a new
     * binary parent never reverses an edge and therefore preserves the input root.
     */
    private static Tree completeTreeRootPreserving(Tree tree, SortedRows rows, int n) {
        boolean[] inTree = new boolean[n];
        TreeNode[] taxonNode = new TreeNode[n];
        TreeNode root = tree.hasPolytomy
            ? deepCopyNodesPolytomy(tree.root, null, taxonNode)
            : deepCopyNodes(tree.root, null, taxonNode);
        for (int i = 0; i < n; i++) if (tree.positionMap[i] != -1) inTree[i] = true;

        for (int x = 0; x < n; x++) {
            if (inTree[x]) continue;
            int anchor = findAnchor(x, inTree, rows);
            TreeNode anchorLeaf = taxonNode[anchor];
            if (anchorLeaf == null) {
                throw new IllegalStateException("Missing completion anchor node for taxon " + anchor);
            }
            TreeNode oldParent = anchorLeaf.parent;
            TreeNode newLeaf = new TreeNode();
            newLeaf.taxonId = x;
            TreeNode joined = new TreeNode();
            setBinaryChildren(joined, anchorLeaf, newLeaf);
            anchorLeaf.parent = joined;
            newLeaf.parent = joined;
            joined.parent = oldParent;
            if (oldParent == null) root = joined;
            else replaceChild(oldParent, anchorLeaf, joined);
            inTree[x] = true;
            taxonNode[x] = newLeaf;
        }

        return tree.hasPolytomy
            ? rebuildTreePolytomy(tree.treeIndex, root, n)
            : rebuildTree(tree.treeIndex, root, n);
    }

    // ── Anchor finding ────────────────────────────────────────────────────────

    /**
     * Find the closest taxon to x that is currently in the tree: the first
     * in-tree taxon in x's nearest-neighbour order.
     *
     * Mirrors AbstractMatrix.getClosestPresentTaxonId() lines 52–67.
     * The original checks (missingId > other || presentBS.get(other)); here we
     * unify both conditions into inTree[candidate] which is true for both
     * originally-present and already-inserted taxa.
     */
    static int findAnchor(int x, boolean[] inTree, SortedRows rows) {
        int n = rows.n;
        for (int rank = 0; rank < n; rank++) {
            int candidate = rows.get(x, rank);
            if (candidate != x && inTree[candidate]) {
                return candidate;
            }
        }
        throw new RuntimeException("No anchor found for taxon " + x
                + " — inTree array may be empty");
    }

    // ── Deep copy ────────────────────────────────────────────────────────────

    /**
     * Recursively deep-copy a TreeNode subtree, also populating taxonNode[].
     *
     * The copies are fresh objects with the same taxonId/rangeStart/rangeEnd
     * but independent parent/left/right pointers — mutations to the copy do
     * not affect the original Tree's nodes (preserving originalTrees' range
     * fields used by PrefixHashArrays and PartitionTable).
     *
     * taxonNode[id] is set to the leaf copy for each leaf encountered.
     *
     * @param src       source node
     * @param parent    parent of the copy (null for root)
     * @param taxonNode map from taxon ID → leaf node; populated for leaves
     * @return root of the copied subtree
     */
    private static TreeNode deepCopyNodes(TreeNode src, TreeNode parent,
                                          TreeNode[] taxonNode) {
        if (src == null) return null;
        TreeNode copy    = new TreeNode();
        copy.taxonId     = src.taxonId;
        copy.isDuplicationNode = src.isDuplicationNode;
        copy.isSpeciationNode = src.isSpeciationNode;
        copy.rangeStart  = src.rangeStart;
        copy.rangeEnd    = src.rangeEnd;
        copy.parent      = parent;
        copy.left        = deepCopyNodes(src.left,  copy, taxonNode);
        copy.right       = deepCopyNodes(src.right, copy, taxonNode);
        if (copy.isLeaf()) {
            taxonNode[copy.taxonId] = copy;
        }
        return copy;
    }

    private static TreeNode deepCopyNodesPolytomy(TreeNode src, TreeNode parent,
                                                   TreeNode[] taxonNode) {
        if (src == null) return null;
        TreeNode copy = new TreeNode();
        copy.taxonId = src.taxonId;
        copy.isDuplicationNode = src.isDuplicationNode;
        copy.isSpeciationNode = src.isSpeciationNode;
        copy.rangeStart = src.rangeStart;
        copy.rangeEnd = src.rangeEnd;
        copy.parent = parent;
        if (src.isLeaf()) {
            taxonNode[copy.taxonId] = copy;
        } else if (src.isPolytomous()) {
            TreeNode[] children = new TreeNode[src.children.length];
            for (int i = 0; i < children.length; i++) {
                children[i] = deepCopyNodesPolytomy(src.children[i], copy, taxonNode);
            }
            setChildren(copy, children);
        } else {
            copy.left = deepCopyNodesPolytomy(src.left, copy, taxonNode);
            copy.right = deepCopyNodesPolytomy(src.right, copy, taxonNode);
        }
        return copy;
    }

    // ── Rebuild Tree ──────────────────────────────────────────────────────────

    /** Reconstruct a Tree object from the mutated TreeNode structure. */
    private static Tree rebuildTree(int treeIndex, TreeNode root, int n) {
        int[] postorderArray = new int[n];
        int[] counter        = {0};
        assignRangesAndFill(root, postorderArray, counter);
        int leafCount = counter[0];
        postorderArray = Arrays.copyOf(postorderArray, leafCount);

        int[] positionMap = new int[n];
        Arrays.fill(positionMap, -1);
        for (int j = 0; j < leafCount; j++) positionMap[postorderArray[j]] = j;

        return new Tree(treeIndex, root, postorderArray, positionMap, leafCount, n);
    }

    private static Tree rebuildTreePolytomy(int treeIndex, TreeNode root, int n) {
        int[] postorderArray = new int[n];
        int[] counter = {0};
        assignRangesAndFillPolytomy(root, postorderArray, counter);
        int leafCount = counter[0];
        postorderArray = Arrays.copyOf(postorderArray, leafCount);
        int[] positionMap = new int[n];
        Arrays.fill(positionMap, -1);
        for (int j = 0; j < leafCount; j++) positionMap[postorderArray[j]] = j;
        return new Tree(treeIndex, root, postorderArray, positionMap, leafCount, n, true);
    }

    private static void assignRangesAndFill(TreeNode node, int[] arr, int[] counter) {
        if (node.isLeaf()) {
            node.rangeStart  = counter[0];
            node.rangeEnd    = counter[0] + 1;
            arr[counter[0]]  = node.taxonId;
            counter[0]++;
            return;
        }
        assignRangesAndFill(node.left,  arr, counter);
        assignRangesAndFill(node.right, arr, counter);
        node.rangeStart = node.left.rangeStart;
        node.rangeEnd   = node.right.rangeEnd;
    }

    private static void assignRangesAndFillPolytomy(TreeNode node, int[] arr, int[] counter) {
        if (node.isLeaf()) {
            node.rangeStart = counter[0];
            node.rangeEnd = counter[0] + 1;
            arr[counter[0]++] = node.taxonId;
            return;
        }
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) {
                assignRangesAndFillPolytomy(child, arr, counter);
            }
        } else {
            assignRangesAndFillPolytomy(node.left, arr, counter);
            assignRangesAndFillPolytomy(node.right, arr, counter);
        }
        node.rangeStart = node.left.rangeStart;
        node.rangeEnd = node.right.rangeEnd;
    }

    // ── Child-list helpers ───────────────────────────────────────────────────

    private static void setChildren(TreeNode node, TreeNode[] children) {
        if (children.length < 2) {
            throw new IllegalArgumentException("Internal node must have at least two children");
        }
        node.left = children[0];
        node.right = children[children.length - 1];
        node.children = children.length > 2 ? children : null;
    }

    private static void setBinaryChildren(TreeNode node, TreeNode left, TreeNode right) {
        node.left = left;
        node.right = right;
        node.children = null;
    }

    private static void replaceChild(TreeNode node, TreeNode remove, TreeNode replacement) {
        if (!node.isPolytomous()) {
            if (node.left == remove) node.left = replacement;
            else if (node.right == remove) node.right = replacement;
            else throw new IllegalStateException(
                "Completion child-list mutation lost the old root");
            return;
        }
        TreeNode[] children = node.children;
        for (int i = 0; i < children.length; i++) {
            if (children[i] == remove) {
                children[i] = replacement;
                setChildren(node, children);
                return;
            }
        }
        throw new IllegalStateException("Completion child-list mutation lost the old root");
    }
}
