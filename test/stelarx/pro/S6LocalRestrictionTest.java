package stelarx.pro;

import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;

import java.lang.reflect.Method;
import java.util.Map;

/** Temporary invariant check for one-pass, disjoint-group tree restriction. */
public final class S6LocalRestrictionTest {
    public static void main(String[] args) throws Exception {
        TaxonRegistry taxa = new TaxonRegistry();
        for (char c = 'a'; c <= 'h'; c++) taxa.register(Character.toString(c));
        taxa.lock();
        Tree tree = TreeParser.parseSpeciesTreeNewick(
            "(((a,b),(c,d)),((e,f),(g,h)));", taxa);
        int[] group = {0, 1, 0, 1, 0, 1, 0, 1};
        int[] local = {0, 0, 1, 1, 2, 2, 3, 3};
        Method method = S6S2BBGuideData.class.getDeclaredMethod(
            "distribute", TreeNode.class, int[].class, int[].class);
        method.setAccessible(true);
        @SuppressWarnings("unchecked")
        Map<Integer, TreeNode> roots =
            (Map<Integer, TreeNode>) method.invoke(null, tree.root, group, local);
        if (roots.size() != 2) throw new AssertionError("Expected two groups");
        for (int g = 0; g < 2; g++) {
            String actual = canonical(roots.get(g));
            if (!actual.equals("((0,1),(2,3))"))
                throw new AssertionError("group " + g + ": " + actual);
        }
        System.out.println("S6 local one-pass restriction: OK");
    }

    private static String canonical(TreeNode node) {
        if (node.isLeaf()) return Integer.toString(node.taxonId);
        if (node.left.parent != node || node.right.parent != node)
            throw new AssertionError("Broken parent link");
        String left = canonical(node.left), right = canonical(node.right);
        return left.compareTo(right) <= 0
            ? "(" + left + "," + right + ")"
            : "(" + right + "," + left + ")";
    }
}
