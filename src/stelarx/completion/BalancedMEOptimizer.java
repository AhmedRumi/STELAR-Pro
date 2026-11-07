package stelarx.completion;

import stelarx.Logging;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.List;

/**
 * Balanced minimum evolution (BME) topology search by NNI hill climbing, seeded from
 * a Neighbor-Joining tree. This is FastME's core loop.
 *
 * <h3>Why this and not BIONJ</h3>
 * BIONJ was measured to produce essentially the same tree as NJ (mean Jaccard 0.969
 * over rooted clades, identical on one condition), because it only changes the
 * <em>reduction formula</em> during agglomeration - a local change that rarely flips a
 * join. BME instead keeps NJ's tree as a starting point and searches topology space
 * under a different objective, so it can move clades that no agglomerative pass would.
 *
 * <h3>The objective</h3>
 * Pauplin's formula gives the balanced tree length of an unrooted binary topology in
 * closed form from the distance matrix alone, with no branch-length fitting:
 *
 * <pre>   L(T) = sum over i&lt;j of  2^(1 - p_ij) * d_ij   </pre>
 *
 * where {@code p_ij} is the number of edges on the path from leaf i to leaf j. On an
 * additive matrix this returns exactly the sum of branch lengths, which is the
 * property the implementation is checked against in
 * {@link #pauplinLengthForTest(double[], int, int[][])}.
 *
 * <h3>Search</h3>
 * First-improvement passes over all {@code n-3} internal edges, two NNI moves each,
 * repeating until a pass makes no change or {@code maxPasses} is reached. Each
 * evaluation recomputes all topological distances, so a pass costs O(n^3); at n = 201
 * that is a few seconds at most against the 10^2-10^3 s the downstream DP takes.
 *
 * <h3>Rooting</h3>
 * NNI works on the unrooted topology, but the DP needs a rooted tree. The optimized
 * topology is rooted on the edge that splits the taxa most evenly. That is a
 * deterministic choice, not a principled one - rooting changes which clusters are
 * donated, so it is a free parameter worth revisiting.
 */
public final class BalancedMEOptimizer {

    private BalancedMEOptimizer() {}

    /** Run NJ, then improve its topology under BME. */
    public static Tree build(SimilarityMatrix sim, int treeIndex,
                             NeighborJoiningClusterer.Transform transform, int maxPasses) {
        Tree seed = NeighborJoiningClusterer.build(
            sim, treeIndex, NeighborJoiningClusterer.Variant.NJ, transform);
        if (seed == null) return null;
        int n = sim.n;
        if (n < 4) return seed;
        double[] d = NeighborJoiningClusterer.distanceMatrix(sim, transform);
        return optimize(d, n, seed, treeIndex, transform, maxPasses);
    }

    static Tree optimize(double[] d, int n, Tree seed, int treeIndex,
                         NeighborJoiningClusterer.Transform transform, int maxPasses) {
        int[][] adj = fromRooted(seed.root, n);
        double before = pauplinLength(d, n, adj);

        int[] internalEdgesU = new int[2 * n], internalEdgesV = new int[2 * n];
        int applied = 0, passes = 0;
        double current = before;

        for (; passes < maxPasses; passes++) {
            int edgeCount = collectInternalEdges(adj, n, internalEdgesU, internalEdgesV);
            boolean improvedThisPass = false;
            for (int e = 0; e < edgeCount; e++) {
                int u = internalEdgesU[e], v = internalEdgesV[e];
                // An earlier move in this pass may have destroyed this edge: NNI keeps
                // its own u-v edge but can move a neighbour from u over to v.
                if (!contains(adj[u], v)) continue;
                // u's two neighbours other than v, v's two other than u
                int[] uo = others(adj[u], v), vo = others(adj[v], u);
                for (int pick = 0; pick < 2; pick++) {
                    int x = uo[1];               // always move u's second neighbour
                    int y = vo[pick];            // against each of v's two in turn
                    swap(adj, u, x, v, y);
                    double after = pauplinLength(d, n, adj);
                    if (after < current - 1e-12) {
                        current = after;
                        applied++;
                        improvedThisPass = true;
                        break;                   // first improvement: keep and move on
                    }
                    swap(adj, u, y, v, x);       // revert
                }
            }
            if (!improvedThisPass) break;
        }

        Logging.info("BME/%s: NNI search applied %d move(s) over %d pass(es), "
                + "balanced tree length %.6f -> %.6f (%.2f%% shorter)",
            transform, applied, passes + 1, before, current,
            before == 0.0 ? 0.0 : 100.0 * (before - current) / before);

        return rootBalanced(adj, n, treeIndex);
    }

    // ── Pauplin balanced tree length ─────────────────────────────────────────

    /** L(T) = sum_{i<j} 2^(1-p_ij) d_ij, with p_ij the topological path length. */
    private static double pauplinLength(double[] d, int n, int[][] adj) {
        int nodes = adj.length;
        int[] dist = new int[nodes];
        ArrayDeque<Integer> queue = new ArrayDeque<>();
        double total = 0.0;
        for (int i = 0; i < n; i++) {
            Arrays.fill(dist, -1);
            dist[i] = 0;
            queue.add(i);
            while (!queue.isEmpty()) {
                int cur = queue.poll();
                for (int nb : adj[cur]) {
                    if (dist[nb] < 0) { dist[nb] = dist[cur] + 1; queue.add(nb); }
                }
            }
            for (int j = i + 1; j < n; j++) {
                total += Math.scalb(d[i * n + j], 1 - dist[j]);   // d * 2^(1-p)
            }
        }
        return total;
    }

    /** Exposed for the self-test: additive input must return the branch-length sum. */
    static double pauplinLengthForTest(double[] d, int n, int[][] adj) {
        return pauplinLength(d, n, adj);
    }

    /** Exposed for the self-test: balanced length of a rooted tree's topology. */
    static double lengthOfRootedForTest(double[] d, int n, Tree t) {
        return pauplinLength(d, n, fromRooted(t.root, n));
    }

    // ── Unrooted topology handling ───────────────────────────────────────────

    /**
     * Convert a rooted binary tree to unrooted adjacency. Leaf node ids equal their
     * taxon ids; internal nodes take n..2n-3. The rooted root is a degree-2 node that
     * has no place in an unrooted binary topology, so it is dropped and its two
     * children joined directly.
     */
    private static int[][] fromRooted(TreeNode root, int n) {
        IdentityHashMap<TreeNode, Integer> id = new IdentityHashMap<>();
        List<int[]> edges = new ArrayList<>(2 * n);
        int[] nextInternal = {n};
        assignIds(root, root, id, nextInternal);
        collectEdges(root, root, id, edges);
        // Join the root's two children, replacing the dropped degree-2 root.
        edges.add(new int[]{id.get(root.left), id.get(root.right)});

        int nodes = 2 * n - 2;
        int[] deg = new int[nodes];
        for (int[] e : edges) { deg[e[0]]++; deg[e[1]]++; }
        int[][] adj = new int[nodes][];
        for (int v = 0; v < nodes; v++) adj[v] = new int[deg[v]];
        int[] fill = new int[nodes];
        for (int[] e : edges) {
            adj[e[0]][fill[e[0]]++] = e[1];
            adj[e[1]][fill[e[1]]++] = e[0];
        }
        return adj;
    }

    private static void assignIds(TreeNode node, TreeNode root,
                                  IdentityHashMap<TreeNode, Integer> id, int[] next) {
        if (node.isLeaf()) { id.put(node, node.taxonId); return; }
        if (node != root) id.put(node, next[0]++);
        assignIds(node.left, root, id, next);
        assignIds(node.right, root, id, next);
    }

    private static void collectEdges(TreeNode node, TreeNode root,
                                     IdentityHashMap<TreeNode, Integer> id,
                                     List<int[]> edges) {
        if (node.isLeaf()) return;
        if (node != root) {
            edges.add(new int[]{id.get(node), id.get(node.left)});
            edges.add(new int[]{id.get(node), id.get(node.right)});
        }
        collectEdges(node.left, root, id, edges);
        collectEdges(node.right, root, id, edges);
    }

    private static int collectInternalEdges(int[][] adj, int n, int[] us, int[] vs) {
        int count = 0;
        for (int u = n; u < adj.length; u++) {
            for (int v : adj[u]) {
                if (v >= n && v > u) { us[count] = u; vs[count] = v; count++; }
            }
        }
        return count;
    }

    private static boolean contains(int[] arr, int value) {
        for (int a : arr) if (a == value) return true;
        return false;
    }

    private static int[] others(int[] neighbours, int exclude) {
        int[] out = new int[neighbours.length - 1];
        int k = 0;
        for (int nb : neighbours) if (nb != exclude) out[k++] = nb;
        return out;
    }

    /** Detach x from u and y from v, then reattach crosswise. */
    private static void swap(int[][] adj, int u, int x, int v, int y) {
        replace(adj[u], x, y);
        replace(adj[v], y, x);
        replace(adj[x], u, v);
        replace(adj[y], v, u);
    }

    private static void replace(int[] arr, int from, int to) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == from) { arr[i] = to; return; }
        }
        throw new IllegalStateException("NNI adjacency corrupted: " + from + " not found");
    }

    // ── Rooting and Tree assembly ────────────────────────────────────────────

    /** Root on the edge whose removal splits the taxa most evenly. */
    private static Tree rootBalanced(int[][] adj, int n, int treeIndex) {
        int bestU = -1, bestV = -1, bestSkew = Integer.MAX_VALUE;
        int nodes = adj.length;
        for (int u = 0; u < nodes; u++) {
            for (int v : adj[u]) {
                if (v <= u) continue;
                int side = leavesBeyond(adj, n, u, v);   // leaves on v's side
                int skew = Math.abs(n - 2 * side);
                if (skew < bestSkew) { bestSkew = skew; bestU = u; bestV = v; }
            }
        }
        TreeNode root = new TreeNode();
        root.left = buildSubtree(adj, n, bestV, bestU);
        root.right = buildSubtree(adj, n, bestU, bestV);
        root.left.parent = root;
        root.right.parent = root;

        int[] postArr = new int[n];
        int[] cursor = {0};
        assignRanges(root, postArr, cursor);
        int[] posMap = new int[n];
        Arrays.fill(posMap, -1);
        for (int pos = 0; pos < n; pos++) posMap[postArr[pos]] = pos;
        return new Tree(treeIndex, root, postArr, posMap, n, n);
    }

    /** Leaves reachable from {@code node} without going back through {@code from}. */
    private static int leavesBeyond(int[][] adj, int n, int from, int node) {
        int count = node < n ? 1 : 0;
        ArrayDeque<int[]> stack = new ArrayDeque<>();
        stack.push(new int[]{node, from});
        while (!stack.isEmpty()) {
            int[] cur = stack.pop();
            for (int nb : adj[cur[0]]) {
                if (nb == cur[1]) continue;
                if (nb < n) count++;
                else stack.push(new int[]{nb, cur[0]});
            }
        }
        return count;
    }

    /** Build the rooted subtree hanging off {@code node}, entered from {@code from}. */
    private static TreeNode buildSubtree(int[][] adj, int n, int node, int from) {
        if (node < n) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId = node;
            return leaf;
        }
        int[] kids = others(adj[node], from);
        TreeNode inner = new TreeNode();
        inner.left = buildSubtree(adj, n, kids[0], node);
        inner.right = buildSubtree(adj, n, kids[1], node);
        inner.left.parent = inner;
        inner.right.parent = inner;
        return inner;
    }

    private static void assignRanges(TreeNode node, int[] postArr, int[] cursor) {
        if (node.isLeaf()) {
            node.rangeStart = cursor[0];
            postArr[cursor[0]++] = node.taxonId;
            node.rangeEnd = cursor[0];
        } else {
            node.rangeStart = cursor[0];
            assignRanges(node.left, postArr, cursor);
            assignRanges(node.right, postArr, cursor);
            node.rangeEnd = cursor[0];
        }
    }
}
