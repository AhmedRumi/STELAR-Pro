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
import stelarx.util.Threading;
import stelarx.weight.WeightTable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

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

    public static List<Tree> generate(S6ScoreSource support, List<Tree> originalGenes,
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

            try (S6ScoreSource scorer = support.fork()) {
            Node root;
            if (n >= 100) {
                int size = Math.max(4, Math.min(n,
                    (int) (Math.sqrt(n) * Math.log(n) / Math.log(2.0) / 4.0)));
                List<Integer> selected = order.subList(0, size);
                if (previousTree == null) {
                    // A-Pro scores subset guides on the induced gene families.
                    // Restriction preserves every triplet among selected taxa
                    // while making each reusable guide state much smaller.
                    TaxonRegistry guideRegistry = new TaxonRegistry();
                    for (int taxon : selected)
                        guideRegistry.register(registry.getName(taxon));
                    guideRegistry.lock();
                    List<Tree> guideGenes = TreeParser.parseRestrictedGeneTrees(
                        taggedInput, guideRegistry, false).trees();
                    // Each guide has its own RNG: parallel workers cannot perturb
                    // proposal order or tie-breaking by sharing the round RNG.
                    List<List<Integer>> guideOrders = new ArrayList<>(size);
                    long[] guideSeeds = new long[size];
                    for (int guide = 0; guide < size; guide++) {
                        List<Integer> guideOrder = new ArrayList<>(size);
                        for (int taxon = 0; taxon < size; taxon++)
                            guideOrder.add(taxon);
                        guideOrders.add(guideOrder);
                        guideSeeds[guide] = random.nextLong();
                    }
                    String[] guideNewicks = new String[size];
                    // Small induced scorers let every guide borrow a CPU state.
                    try (S6DynamicTripletSupport guideSupport =
                             S6DynamicTripletSupport.build(guideGenes, size,
                                 Config.getInstance().getThreadCount())) {
                        int cpuThreads = Threading.isStarted()
                            ? Threading.getNumThreads() : 1;
                        Logging.info("S6 round %d/%d: %d triplet guides; up to %d "
                            + "concurrent scorer(s) in a %d-thread CPU pool; "
                            + "backbone uses STELAR triplet DP", round + 1,
                            rounds, size, Math.min(size,
                                Math.min(guideSupport.scoringStates(), cpuThreads)),
                            cpuThreads);
                        int[] guideTaxa = new int[size];
                        for (int taxon = 0; taxon < size; taxon++)
                            guideTaxa[taxon] = taxon;
                        guideSupport.prepareGuideFocus(guideTaxa);
                        if (Threading.isStarted() && Threading.getNumThreads() > 1) {
                            Threading.processRangeParallel(size, guide ->
                                guideNewicks[guide] = buildGuide(guideOrders.get(guide),
                                    guideSeeds[guide], guideSupport, guideRegistry, n));
                        } else {
                            for (int guide = 0; guide < size; guide++)
                                guideNewicks[guide] = buildGuide(guideOrders.get(guide),
                                    guideSeeds[guide], guideSupport, guideRegistry, n);
                        }
                    }
                    List<String> guides = java.util.Arrays.asList(guideNewicks);
                    root = inferBackbone(guideGenes, guideRegistry, guides, registry);
                } else {
                    root = inducedPrevious(previousTree, selected, registry);
                }
                root = placeByBackboneGroups(root, order.subList(size, n),
                    scorer, random);
            } else if (previousTree == null) {
                root = fromOrder(order, scorer, random);
            } else {
                int size = Math.min(n, Math.max(4, n / 4));
                root = inducedPrevious(previousTree, order.subList(0, size), registry);
                for (int i = size; i < n; i++) {
                    Edge edge = locate(root, order.get(i), scorer, -1, null, random);
                    root = insert(root, edge.node, order.get(i), -1);
                }
            }
            refineNni(root, scorer, n, true);
            String treeNewick = newick(root, registry) + ";";
            Tree parsed = TreeParser.parseSpeciesTreeNewick(treeNewick, registry);
            markSpeciation(parsed.root);
            proposals.add(new Tree(firstIndex + round, parsed.root,
                parsed.postorderArray, parsed.positionMap, parsed.leafCount, n));
            // The outer STELAR DP scores these proposals against original
            // genes; rescoring every greedy edge here is diagnostic only.
            Logging.info("S6 triplet-greedy proposal %d/%d constructed",
                round + 1, rounds);
            }
        }
        return proposals;
    }

    private static String buildGuide(List<Integer> order, long seed,
                                     S6ScoreSource support,
                                     TaxonRegistry registry, int fullTaxa) {
        try (S6ScoreSource scorer = support.fork()) {
            Random random = new Random(seed);
            Collections.shuffle(order, random);
            Node tree = fromOrder(order, scorer, random);
            // Guide workers must not submit nested jobs to the same fixed pool.
            refineNni(tree, scorer, fullTaxa, false);
            return newick(tree, registry) + ";";
        }
    }

    /** One focused slice supplies every triplet needed to place taxon x. */
    private static Edge locate(Node root, int x, S6ScoreSource support,
                               int group, Node originalChild, Random random) {
        placementCounts(root, x, support.focus(x));
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

    private static void placementCounts(Node node, int x, S6ScoreSource.Focus support) {
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
        // Only ancestors of the inserted edge change their taxon sets.
        if (edge.taxa != null) {
            joint.taxa = ArraysCopy.join(edge.taxa, added.taxa);
            for (Node up = parent; up != null; up = up.parent)
                up.taxa = ArraysCopy.join(up.left.taxa, up.right.taxa);
        }
        return root;
    }

    /** Test hook: compare every insertion delta against full fixed-tree rescoring. */
    static void verifyInsertionDeltas(S6ScoreSource support) {
        if (support.taxa() < 4) throw new IllegalArgumentException("Need four taxa");
        Node root = new Node(new Node(new Node(0), new Node(1)), new Node(2));
        refresh(root);
        placementCounts(root, 3, support.focus(3));
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

    private static Node fromOrder(List<Integer> order, S6ScoreSource support,
                                  Random random) {
        if (order.size() < 2) throw new IllegalArgumentException("Need two guide taxa");
        Node root = new Node(new Node(order.get(0)), new Node(order.get(1)));
        refresh(root);
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
                                               S6ScoreSource support, Random random) {
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
        refineNni(root, support, support.taxa(), true);
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
    private static void refineNni(Node root, S6ScoreSource support, int n,
                                  boolean allowParallel) {
        int limit = 20 + (int) (2 * Math.sqrt(n) * Math.log(n) / Math.log(2.0));
        for (int move = 0; move < limit; move++) {
            refresh(root);
            List<Node> nodes = new ArrayList<>();
            collectPostorder(root, nodes);
            int index = -1, choice = 0;
            if (allowParallel && n >= 100 && Threading.isStarted()
                    && Threading.getNumThreads() > 1) {
                // Score disjoint candidate NNIs on CPU workers, but apply only
                // the first improving move in postorder, as in serial S6.
                int[] choices = new int[nodes.size()];
                AtomicInteger earliest = new AtomicInteger(nodes.size());
                Threading.processRangeParallel(nodes.size(), i -> {
                    if (i >= earliest.get()) return;
                    try (S6ScoreSource worker = support.fork()) {
                        int candidate = nniChoice(nodes.get(i), worker);
                        choices[i] = candidate;
                        if (candidate != 0) earliest.accumulateAndGet(i, Math::min);
                    }
                });
                index = earliest.get();
                if (index < nodes.size()) choice = choices[index];
            } else {
                for (int i = 0; i < nodes.size(); i++) {
                    choice = nniChoice(nodes.get(i), support);
                    if (choice != 0) { index = i; break; }
                }
            }
            if (index < 0 || index >= nodes.size()) break;
            applyNni(nodes.get(index), choice);
        }
        refresh(root);
    }

    /** Encode left/right child and which of its two arms is moved. */
    private static int nniChoice(Node parent, S6ScoreSource support) {
        if (parent.leaf()) return 0;
        for (int side = 0; side < 2; side++) {
            Node child = side == 0 ? parent.left : parent.right;
            if (child.leaf()) continue;
            Node sibling = side == 0 ? parent.right : parent.left;
            long[] scores = support.groupScores(
                child.left.taxa, child.right.taxa, sibling.taxa);
            long old = scores[0], swapB = scores[1], swapA = scores[2];
            if (swapB > old || swapA > old)
                return 1 + 2 * side + (swapB >= swapA ? 0 : 1);
        }
        return 0;
    }

    private static void applyNni(Node parent, int choice) {
        boolean first = choice <= 2;
        Node child = first ? parent.left : parent.right;
        Node sibling = first ? parent.right : parent.left;
        Node a = child.left, b = child.right;
        Node moved = (choice & 1) == 1 ? b : a;
        Node kept = (choice & 1) == 1 ? a : b;
        child.left = kept;
        child.right = sibling;
        kept.parent = sibling.parent = child;
        if (first) parent.right = moved;
        else parent.left = moved;
        moved.parent = parent;
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

    private static long score(Node node, S6ScoreSource support) {
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

    /** Match S5: guides supply subset candidates; restricted genes supply weights. */
    private static Node inferBackbone(List<Tree> genes, TaxonRegistry subset,
                                      List<String> guideNewicks,
                                      TaxonRegistry registry) throws IOException {
        List<Tree> guides = new ArrayList<>(guideNewicks.size());
        for (String newick : guideNewicks) {
            Tree parsed = TreeParser.parseSpeciesTreeNewick(newick, subset);
            markSpeciation(parsed.root);
            guides.add(new Tree(guides.size(), parsed.root,
                parsed.postorderArray, parsed.positionMap,
                parsed.leafCount, subset.size()));
        }
        TaxonHasher hasher = new TaxonHasher(subset.size(),
            Config.getInstance().getNumHashSeeds(), Config.getInstance().getBaseSeed());
        PrefixHashArrays candidatePrefix = new PrefixHashArrays(guides, hasher);
        PrefixHashArrays genePrefix = new PrefixHashArrays(genes, hasher);
        UniqueTaxonSubtreeHashes candidateHashes =
            new UniqueTaxonSubtreeHashes(guides, hasher);
        UniqueTaxonSubtreeHashes geneHashes =
            new UniqueTaxonSubtreeHashes(genes, hasher);
        ClusterTable clusters = new ClusterTable(guides, candidatePrefix,
            subset.size(), candidateHashes);
        PartitionTable partitions = new PartitionTable(genes, genePrefix, geneHashes);
        DPTable dp = new DPTable(guides, candidatePrefix, clusters,
            candidateHashes, Integer.MAX_VALUE);
        candidateHashes.release();
        geneHashes.release();
        // Keep subset-backbone recombination consistent with the outer DP mode.
        if (Config.getInstance().getSearchMode() == Config.SearchMode.FULL)
            dp.addCrossTreeTransitions(clusters,
                Config.getInstance().getComputeMode() == Config.ComputeMode.GPU
                    && GPUDPBuilder.tryLoad());
        WeightTable weights = new WeightTable(dp, partitions, clusters,
            guides, genes);
        String inferred = new Inference().run(dp, weights, clusters,
            guides, subset, hasher);
        Node result = fromSpeciesTree(
            TreeParser.parseSpeciesTreeNewick(inferred, subset).root,
            subset, registry);
        refresh(result);
        return result;
    }
}
