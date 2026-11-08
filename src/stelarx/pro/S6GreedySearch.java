package stelarx.pro;

import stelarx.Config;
import stelarx.Logging;
import stelarx.cluster.ClusterTable;
import stelarx.completion.SimilarityMatrix;
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
import java.util.PriorityQueue;
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
    private record GroupPair(double score, int a, int b) {}

    public static List<Tree> generate(S6ScoreSource support, List<Tree> originalGenes,
                                      String taggedInput, TaxonRegistry registry,
                                      int firstIndex, long seed, int rounds,
                                      String previousTree) throws IOException {
        return generate(support, originalGenes, taggedInput, registry,
            firstIndex, seed, rounds, previousTree, null);
    }

    /** Shared guide data supplies S4-style backbones and UPGMA partitions. */
    public static List<Tree> generate(S6ScoreSource support, List<Tree> originalGenes,
                                      String taggedInput, TaxonRegistry registry,
                                      int firstIndex, long seed, int rounds,
                                      String previousTree,
                                      S6S2BBGuideData guideData) throws IOException {
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
                if (previousTree == null && guideData != null) {
                    TaxonRegistry subset = new TaxonRegistry();
                    for (int taxon : selected)
                        subset.register(registry.getName(taxon));
                    subset.lock();
                    List<Tree> genes = TreeParser.parseRestrictedGeneTrees(
                        taggedInput, subset, false).trees();
                    root = inferBackboneTrees(genes, subset,
                        guideData.restrict(subset), registry);
                    Logging.info("S6-S2BB round %d/%d: STELAR triplet-DP backbone "
                        + "from restricted S4 candidates", round + 1, rounds);
                } else if (previousTree == null) {
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
                root = guideData == null
                    ? placeByBackboneGroups(root, order.subList(size, n),
                        scorer, random)
                    : placeByBackboneGroupsUpgma(root, order.subList(size, n),
                        scorer, guideData.matrix(), random);
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
        return locate(root, x, support, group, originalChild, random, null);
    }

    private static Edge locate(Node root, int x, S6ScoreSource support,
                               int group, Node originalChild, Random random,
                               int[] backboneTaxa) {
        placementCounts(root, x, backboneTaxa == null
            ? support.focus(x) : support.focusOn(x, backboneTaxa));
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
            refresh(copied); // copy() links nodes but does not populate internal sets.
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

    /** Triplet-scored edge assignment; UPGMA resolves each edge's taxa once. */
    private static Node placeByBackboneGroupsUpgma(Node root,
                                                    List<Integer> remaining,
                                                    S6ScoreSource support,
                                                    SimilarityMatrix matrix,
                                                    Random random) {
        List<Node> edges = new ArrayList<>();
        collectNodes(root, edges);
        IdentityHashMap<Node, Integer> edgeId = new IdentityHashMap<>();
        List<List<Integer>> groups = new ArrayList<>(edges.size());
        for (int i = 0; i < edges.size(); i++) {
            edgeId.put(edges.get(i), i);
            groups.add(new ArrayList<>());
        }
        // The backbone is unchanged while assigning taxa, so each focused
        // triplet placement score needs only backbone-species pairs.
        int[] backboneTaxa = root.taxa;
        for (int taxon : remaining) {
            Node edge = locate(root, taxon, support, -1, null, random,
                backboneTaxa).node;
            groups.get(edgeId.get(edge)).add(taxon);
        }
        int nonempty = 0;
        for (int i = 0; i < edges.size(); i++) {
            if (groups.get(i).isEmpty()) continue;
            nonempty++;
            Node edge = edges.get(i);
            Node parent = edge.parent;
            Node joined = new Node(edge, upgmaGroup(groups.get(i), matrix));
            joined.parent = parent;
            if (parent == null) root = joined;
            else if (parent.left == edge) parent.left = joined;
            else if (parent.right == edge) parent.right = joined;
            else throw new IllegalStateException("S6-S2BB backbone parent changed");
        }
        refresh(root);
        Logging.info("S6-S2BB UPGMA placement: %d taxa in %d nonempty "
            + "backbone-edge groups", remaining.size(), nonempty);
        return root;
    }

    /** Group-sized scratch only; the packed global similarity stays read-only. */
    private static Node upgmaGroup(List<Integer> group, SimilarityMatrix global) {
        int m = group.size();
        if (m == 1) return new Node(group.get(0));
        int capacity = 2 * m - 1;
        double[][] sim = new double[capacity][capacity];
        Node[] roots = new Node[capacity];
        int[] sizes = new int[capacity];
        boolean[] active = new boolean[capacity];
        PriorityQueue<GroupPair> nearest = new PriorityQueue<>((x, y) -> {
            int cmp = Double.compare(y.score(), x.score());
            if (cmp != 0) return cmp;
            cmp = Integer.compare(y.a(), x.a());
            return cmp != 0 ? cmp : Integer.compare(y.b(), x.b());
        });
        for (int a = 0; a < m; a++) {
            roots[a] = new Node(group.get(a));
            sizes[a] = 1;
            active[a] = true;
            for (int b = 0; b < a; b++) {
                double score = global.getSim(group.get(a), group.get(b));
                sim[a][b] = sim[b][a] = score;
                nearest.add(new GroupPair(score, b, a));
            }
        }
        for (int next = m; next < capacity; next++) {
            GroupPair best;
            do {
                best = nearest.poll();
                if (best == null)
                    throw new IllegalStateException("S6-S2BB UPGMA exhausted pairs");
            } while (!active[best.a()] || !active[best.b()]);
            int a = best.a(), b = best.b();
            roots[next] = new Node(roots[a], roots[b]);
            sizes[next] = sizes[a] + sizes[b];
            active[a] = active[b] = false;
            for (int k = 0; k < next; k++) {
                if (!active[k]) continue;
                double score = (sizes[a] * sim[a][k] + sizes[b] * sim[b][k])
                    / sizes[next];
                sim[next][k] = sim[k][next] = score;
                nearest.add(new GroupPair(score, k, next));
            }
            active[next] = true;
        }
        return roots[capacity - 1];
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
        int parallelism = allowParallel && n >= 100 && Threading.isStarted()
            ? Math.min(3, Threading.getNumThreads()) : 1;
        long started = System.nanoTime();
        // An NNI changes only the swapped neighborhood, so choices elsewhere
        // remain valid even though their position in the postorder may move.
        IdentityHashMap<Node, Integer> cachedChoices = new IdentityHashMap<>();
        if (n >= 1000)
            Logging.info("S6 triplet NNI: up to %d first-improving moves, "
                + "%d bounded scorer worker(s)", limit, parallelism);
        int applied = 0;
        for (int move = 0; move < limit; move++) {
            refresh(root);
            List<Node> nodes = new ArrayList<>();
            collectPostorder(root, nodes);
            int index = -1, choice = 0;
            if (parallelism > 1) {
                // A full-tree task per node previously queued thousands of
                // jobs against just three colour states. Scan small postorder
                // windows, borrowing each state once per worker; the minimum
                // improving index is exactly the serial first-improvement.
                int window = 24 * parallelism;
                for (int begin = 0; begin < nodes.size(); begin += window) {
                    final int from = begin;
                    final int end = Math.min(nodes.size(), begin + window);
                    final int chunk = (end - from + parallelism - 1) / parallelism;
                    int[] choices = new int[end - from];
                    int earliestKnown = end;
                    for (int i = from; i < end; i++) {
                        Integer cached = cachedChoices.get(nodes.get(i));
                        choices[i - from] = cached == null
                            ? Integer.MIN_VALUE : cached;
                        if (cached != null && cached != 0 && i < earliestKnown)
                            earliestKnown = i;
                    }
                    boolean missingBeforeEarliest = false;
                    for (int i = from; i < earliestKnown; i++) {
                        if (choices[i - from] == Integer.MIN_VALUE) {
                            missingBeforeEarliest = true;
                            break;
                        }
                    }
                    if (!missingBeforeEarliest) {
                        if (earliestKnown < end) {
                            index = earliestKnown;
                            choice = choices[index - from];
                            break;
                        }
                        continue;
                    }
                    AtomicInteger earliest = new AtomicInteger(earliestKnown);
                    Threading.processRangeParallel(parallelism, workerId -> {
                        S6ScoreSource worker = workerId == 0 ? support : support.fork();
                        try {
                            int first = from + workerId * chunk;
                            int last = Math.min(end, first + chunk);
                            for (int i = first; i < last && i < earliest.get(); i++) {
                                if (choices[i - from] != Integer.MIN_VALUE)
                                    continue;
                                int candidate = nniChoice(nodes.get(i), worker);
                                choices[i - from] = candidate;
                                if (candidate != 0)
                                    earliest.accumulateAndGet(i, Math::min);
                            }
                        } finally {
                            if (workerId != 0) worker.close();
                        }
                    });
                    // Workers only read the cache; publish their results after
                    // joining so the IdentityHashMap is never mutated in parallel.
                    for (int i = from; i < end; i++) {
                        int candidate = choices[i - from];
                        if (candidate != Integer.MIN_VALUE)
                            cachedChoices.put(nodes.get(i), candidate);
                    }
                    if (earliest.get() < end) {
                        index = earliest.get();
                        choice = choices[index - from];
                        break;
                    }
                }
            } else {
                for (int i = 0; i < nodes.size(); i++) {
                    Node node = nodes.get(i);
                    Integer cached = cachedChoices.get(node);
                    choice = cached == null ? nniChoice(node, support) : cached;
                    if (cached == null) cachedChoices.put(node, choice);
                    if (choice != 0) { index = i; break; }
                }
            }
            if (index < 0 || index >= nodes.size()) break;
            Node parent = nodes.get(index);
            Node child = choice <= 2 ? parent.left : parent.right;
            Node sibling = choice <= 2 ? parent.right : parent.left;
            Node a = child.left, b = child.right;
            Node grandparent = parent.parent;
            applyNni(parent, choice);
            // The parent keeps the same species set. Only these local NNI
            // choices can change; grandparent observes the parent's new arms.
            cachedChoices.remove(grandparent);
            cachedChoices.remove(parent);
            cachedChoices.remove(child);
            cachedChoices.remove(sibling);
            cachedChoices.remove(a);
            cachedChoices.remove(b);
            applied++;
            if (n >= 1000 && applied % 20 == 0)
                Logging.info("S6 triplet NNI: %d move(s) applied in %d ms",
                    applied, (System.nanoTime() - started) / 1_000_000);
        }
        refresh(root);
        if (n >= 1000)
            Logging.info("S6 triplet NNI: finished %d move(s) in %d ms",
                applied, (System.nanoTime() - started) / 1_000_000);
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
        return inferBackboneTrees(genes, subset, guides, registry);
    }

    private static Node inferBackboneTrees(List<Tree> genes, TaxonRegistry subset,
                                           List<Tree> guides,
                                           TaxonRegistry registry) throws IOException {
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
