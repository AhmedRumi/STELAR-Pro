package stelarx.pro;

import stelarx.Config;
import stelarx.Logging;
import stelarx.cluster.ClusterTable;
import stelarx.dp.DPTable;
import stelarx.dp.Inference;
import stelarx.gpu.GPUDPBuilder;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;
import stelarx.weight.WeightTable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Random;

/** S6's independent, rooted-triplet greedy search; no A-Pro objective is called. */
public final class S6GreedySearch {
    private S6GreedySearch() {}

    private static final class Node {
        final int taxon;
        Node left, right, parent;
        int[] taxa;
        int group = -1; // assigned backbone edge; used only during grouped placement
        long pair, leftDelta, rightDelta, prefix;

        Node(int taxon) { this.taxon = taxon; this.taxa = new int[]{taxon}; }
        Node(Node left, Node right) {
            this.taxon = -1;
            this.left = left;
            this.right = right;
            left.parent = right.parent = this;
        }
        boolean leaf() { return taxon >= 0; }
    }

    private record Edge(Node node, long delta) {}

    public static List<Tree> generate(S6TripletSupport support, List<Tree> originalGenes,
                                      String taggedInput, TaxonRegistry registry,
                                      int firstIndex, long seed, int rounds,
                                      String previousTree) throws IOException {
        int n = registry.size();
        if (support.taxa() != n || rounds < 1)
            throw new IllegalArgumentException("Invalid S6 support or round count");
        List<Tree> proposals = new ArrayList<>(rounds);
        Random random = new Random(seed);
        for (int round = 0; round < rounds; round++) {
            List<Integer> order = new ArrayList<>(n);
            for (int taxon = 0; taxon < n; taxon++) order.add(taxon);
            Collections.shuffle(order, random);

            Node root;
            if (n >= 100) {
                int size = Math.max(4, Math.min(n,
                    (int) (Math.sqrt(n) * Math.log(n) / Math.log(2.0) / 4.0)));
                List<Integer> selected = order.subList(0, size);
                if (previousTree == null) {
                    List<String> guides = new ArrayList<>(size);
                    for (int guide = 0; guide < size; guide++) {
                        List<Integer> guideOrder = new ArrayList<>(selected);
                        Collections.shuffle(guideOrder, random);
                        Node tree = fromOrder(guideOrder, support, random);
                        refineNni(tree, support, n);
                        guides.add(newick(tree, registry) + ";");
                    }
                    root = inferBackbone(taggedInput, selected, guides, registry);
                } else {
                    root = inducedPrevious(previousTree, selected, registry);
                }
                root = placeByBackboneGroups(root, order.subList(size, n),
                    support, random);
            } else if (previousTree == null) {
                root = fromOrder(order, support, random);
            } else {
                int size = Math.min(n, Math.max(4, n / 4));
                root = inducedPrevious(previousTree, order.subList(0, size), registry);
                for (int i = size; i < n; i++) {
                    Edge edge = locate(root, order.get(i), support, -1, null, random);
                    root = insert(root, edge.node, order.get(i), -1);
                }
            }
            refineNni(root, support, n);
            String treeNewick = newick(root, registry) + ";";
            Tree parsed = TreeParser.parseSpeciesTreeNewick(treeNewick, registry);
            markSpeciation(parsed.root);
            proposals.add(new Tree(firstIndex + round, parsed.root,
                parsed.postorderArray, parsed.positionMap, parsed.leafCount, n));
            Logging.info("S6 triplet-greedy proposal %d/%d: score %d",
                round + 1, rounds, score(root, support));
        }
        return proposals;
    }

    /** All insertion edges are scored in O(n²) tensor lookups for one new taxon. */
    private static Edge locate(Node root, int x, S6TripletSupport support,
                               int group, Node originalChild, Random random) {
        refresh(root);
        placementCounts(root, x, support);
        List<Edge> edges = new ArrayList<>();
        collectEdges(root, 0L, edges);
        long best = Long.MIN_VALUE;
        List<Node> ties = new ArrayList<>();
        for (Edge edge : edges) {
            Node node = edge.node;
            if (group >= 0 && node != originalChild && node.group != group) continue;
            if (edge.delta > best) {
                best = edge.delta;
                ties.clear();
                ties.add(node);
            } else if (edge.delta == best) {
                ties.add(node);
            }
        }
        if (ties.isEmpty()) throw new IllegalStateException("No S6 placement edge");
        return new Edge(ties.get(random.nextInt(ties.size())), best);
    }

    private static void placementCounts(Node node, int x, S6TripletSupport support) {
        if (node.leaf()) { node.pair = 0L; return; }
        placementCounts(node.left, x, support);
        placementCounts(node.right, x, support);
        long across = 0L, left = 0L, right = 0L;
        for (int a : node.left.taxa) {
            for (int b : node.right.taxa) {
                across += support.support(a, b, x);
                left += support.support(x, a, b);
                right += support.support(x, b, a);
            }
        }
        node.pair = node.left.pair + node.right.pair + across;
        node.leftDelta = left;
        node.rightDelta = right;
    }

    private static void collectEdges(Node node, long prefix, List<Edge> edges) {
        node.prefix = prefix;
        edges.add(new Edge(node, prefix + node.pair));
        if (node.leaf()) return;
        collectEdges(node.left, prefix + node.leftDelta + node.right.pair, edges);
        collectEdges(node.right, prefix + node.rightDelta + node.left.pair, edges);
    }

    private static Node insert(Node root, Node edge, int x, int group) {
        Node parent = edge.parent;
        Node added = new Node(x);
        Node joint = new Node(edge, added);
        joint.group = added.group = group;
        joint.parent = parent;
        if (parent == null) root = joint;
        else if (parent.left == edge) parent.left = joint;
        else parent.right = joint;
        return root;
    }

    /** Test hook: compare every insertion delta against full fixed-tree rescoring. */
    static void verifyInsertionDeltas(S6TripletSupport support) {
        if (support.taxa() < 4) throw new IllegalArgumentException("Need four taxa");
        Node root = new Node(new Node(new Node(0), new Node(1)), new Node(2));
        refresh(root);
        placementCounts(root, 3, support);
        List<Edge> edges = new ArrayList<>();
        collectEdges(root, 0L, edges);
        long base = score(root, support);
        for (Edge edge : edges) {
            IdentityHashMap<Node, Node> copies = new IdentityHashMap<>();
            Node copied = copy(root, copies);
            copied = insert(copied, copies.get(edge.node), 3, -1);
            refresh(copied);
            long observed = score(copied, support) - base;
            if (observed != edge.delta)
                throw new AssertionError("S6 insertion delta mismatch: expected "
                    + observed + ", predicted " + edge.delta);
        }
    }

    private static Node copy(Node node, IdentityHashMap<Node, Node> copies) {
        Node result = node.leaf()
            ? new Node(node.taxon)
            : new Node(copy(node.left, copies), copy(node.right, copies));
        copies.put(node, result);
        return result;
    }

    private static Node fromOrder(List<Integer> order, S6TripletSupport support,
                                  Random random) {
        if (order.size() < 2) throw new IllegalArgumentException("Need two guide taxa");
        Node root = new Node(new Node(order.get(0)), new Node(order.get(1)));
        for (int i = 2; i < order.size(); i++) {
            int x = order.get(i);
            root = insert(root, locate(root, x, support, -1, null, random).node,
                x, -1);
        }
        refresh(root);
        return root;
    }

    /** Mirror A-Pro's backbone-edge partitioning, followed by global cleanup. */
    private static Node placeByBackboneGroups(Node root, List<Integer> remaining,
                                               S6TripletSupport support, Random random) {
        List<Node> backboneEdges = new ArrayList<>();
        collectNodes(root, backboneEdges);
        List<List<Integer>> groups = new ArrayList<>(backboneEdges.size());
        for (int i = 0; i < backboneEdges.size(); i++) groups.add(new ArrayList<>());
        for (int x : remaining) {
            Node edge = locate(root, x, support, -1, null, random).node;
            groups.get(backboneEdges.indexOf(edge)).add(x);
        }
        List<Integer> groupOrder = new ArrayList<>(backboneEdges.size());
        for (int i = 0; i < backboneEdges.size(); i++) groupOrder.add(i);
        Collections.shuffle(groupOrder, random);
        List<Integer> deferred = new ArrayList<>();
        for (int group : groupOrder) {
            Collections.shuffle(groups.get(group), random);
            Node originalChild = backboneEdges.get(group);
            for (int x : groups.get(group)) {
                Edge global = locate(root, x, support, -1, null, random);
                if (global.node != originalChild && global.node.group != group) {
                    deferred.add(x);
                    continue;
                }
                Edge local = locate(root, x, support, group, originalChild, random);
                root = insert(root, local.node, x, group);
            }
        }
        refineNni(root, support, support.taxa());
        Collections.shuffle(deferred, random);
        for (int x : deferred)
            root = insert(root, locate(root, x, support, -1, null, random).node,
                x, -1);
        Logging.info("S6 two-step placement: %d grouped, %d deferred",
            remaining.size() - deferred.size(), deferred.size());
        return root;
    }

    private static void collectNodes(Node node, List<Node> result) {
        result.add(node);
        if (!node.leaf()) {
            collectNodes(node.left, result);
            collectNodes(node.right, result);
        }
    }

    /** A rooted NNI changes exactly the triplets with one taxon in each arm. */
    private static void refineNni(Node root, S6TripletSupport support, int n) {
        int limit = 20 + (int) (2 * Math.sqrt(n) * Math.log(n) / Math.log(2.0));
        for (int move = 0; move < limit; move++) {
            refresh(root);
            List<Node> nodes = new ArrayList<>();
            collectPostorder(root, nodes);
            boolean changed = false;
            for (Node parent : nodes) {
                if (parent.leaf()) continue;
                for (boolean first : new boolean[]{true, false}) {
                    Node child = first ? parent.left : parent.right;
                    if (child.leaf()) continue;
                    Node sibling = first ? parent.right : parent.left;
                    long old = 0, swapB = 0, swapA = 0;
                    for (int a : child.left.taxa)
                        for (int b : child.right.taxa)
                            for (int c : sibling.taxa) {
                                old += support.support(a, b, c);
                                swapB += support.support(a, c, b);
                                swapA += support.support(b, c, a);
                            }
                    if (swapB > old || swapA > old) {
                        Node a = child.left, b = child.right;
                        Node moved = swapB >= swapA ? b : a;
                        Node kept = swapB >= swapA ? a : b;
                        child.left = kept;
                        child.right = sibling;
                        kept.parent = sibling.parent = child;
                        if (first) parent.right = moved;
                        else parent.left = moved;
                        moved.parent = parent;
                        changed = true;
                        break;
                    }
                }
                if (changed) break;
            }
            if (!changed) break;
        }
        refresh(root);
    }

    private static void collectPostorder(Node node, List<Node> result) {
        if (!node.leaf()) {
            collectPostorder(node.left, result);
            collectPostorder(node.right, result);
        }
        result.add(node);
    }

    private static int[] refresh(Node node) {
        if (node.leaf()) return node.taxa;
        int[] a = refresh(node.left), b = refresh(node.right);
        int[] joined = ArraysCopy.join(a, b);
        node.taxa = joined;
        return joined;
    }

    private static final class ArraysCopy {
        static int[] join(int[] a, int[] b) {
            int[] result = new int[a.length + b.length];
            System.arraycopy(a, 0, result, 0, a.length);
            System.arraycopy(b, 0, result, a.length, b.length);
            return result;
        }
    }

    private static long score(Node node, S6TripletSupport support) {
        if (node.leaf()) return 0;
        return score(node.left, support) + score(node.right, support)
            + support.splitWeight(node.left.taxa, node.right.taxa);
    }

    private static String newick(Node node, TaxonRegistry registry) {
        if (node.leaf()) return registry.getName(node.taxon);
        return "(" + newick(node.left, registry) + ","
            + newick(node.right, registry) + ")";
    }

    private static void markSpeciation(TreeNode node) {
        if (node.isLeaf()) return;
        node.isSpeciationNode = true;
        node.isDuplicationNode = false;
        markSpeciation(node.left);
        markSpeciation(node.right);
    }

    private static Node fromSpeciesTree(TreeNode node, TaxonRegistry from,
                                        TaxonRegistry to) {
        if (node.isLeaf()) return new Node(to.getId(from.getName(node.taxonId)));
        return new Node(fromSpeciesTree(node.left, from, to),
            fromSpeciesTree(node.right, from, to));
    }

    private static Node inducedPrevious(String previous, List<Integer> selected,
                                        TaxonRegistry registry) {
        boolean[] keep = new boolean[registry.size()];
        for (int taxon : selected) keep[taxon] = true;
        Tree tree = TreeParser.parseSpeciesTreeNewick(previous, registry);
        Node result = prune(tree.root, keep);
        if (result == null) throw new IllegalStateException("Empty S6 previous backbone");
        refresh(result);
        return result;
    }

    private static Node prune(TreeNode node, boolean[] keep) {
        if (node.isLeaf()) return keep[node.taxonId] ? new Node(node.taxonId) : null;
        Node left = prune(node.left, keep), right = prune(node.right, keep);
        if (left == null) return right;
        if (right == null) return left;
        return new Node(left, right);
    }

    /** Subset DP uses original restricted genes for weights and triplet guides for X. */
    private static Node inferBackbone(String taggedInput, List<Integer> selected,
                                      List<String> guideNewicks,
                                      TaxonRegistry registry) throws IOException {
        TaxonRegistry subset = new TaxonRegistry();
        for (int taxon : selected) subset.register(registry.getName(taxon));
        subset.lock();
        List<Tree> genes = TreeParser.parseRestrictedGeneTrees(
            taggedInput, subset, false).trees();
        List<Tree> candidates = new ArrayList<>(genes);
        for (String newick : guideNewicks) {
            Tree parsed = TreeParser.parseSpeciesTreeNewick(newick, subset);
            markSpeciation(parsed.root);
            candidates.add(new Tree(candidates.size(), parsed.root,
                parsed.postorderArray, parsed.positionMap,
                parsed.leafCount, subset.size()));
        }
        TaxonHasher hasher = new TaxonHasher(subset.size(),
            Config.getInstance().getNumHashSeeds(), Config.getInstance().getBaseSeed());
        PrefixHashArrays candidatePrefix = new PrefixHashArrays(candidates, hasher);
        PrefixHashArrays genePrefix = new PrefixHashArrays(genes, hasher);
        UniqueTaxonSubtreeHashes candidateHashes =
            new UniqueTaxonSubtreeHashes(candidates, hasher);
        UniqueTaxonSubtreeHashes geneHashes =
            new UniqueTaxonSubtreeHashes(genes, hasher);
        ClusterTable clusters = new ClusterTable(candidates, candidatePrefix,
            subset.size(), candidateHashes);
        PartitionTable partitions = new PartitionTable(genes, genePrefix, geneHashes);
        DPTable dp = new DPTable(candidates, candidatePrefix, clusters,
            candidateHashes, Integer.MAX_VALUE);
        candidateHashes.release();
        geneHashes.release();
        dp.addCrossTreeTransitions(clusters,
            Config.getInstance().getComputeMode() == Config.ComputeMode.GPU
                && GPUDPBuilder.tryLoad());
        WeightTable weights = new WeightTable(dp, partitions, clusters,
            candidates, genes);
        String inferred = new Inference().run(dp, weights, clusters,
            candidates, subset, hasher);
        Node result = fromSpeciesTree(
            TreeParser.parseSpeciesTreeNewick(inferred, subset).root,
            subset, registry);
        refresh(result);
        return result;
    }
}
