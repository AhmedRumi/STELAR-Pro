package stelarx.pro;

import stelarx.cluster.ClusterHash;
import stelarx.cluster.ClusterTable;
import stelarx.dp.BipartitionSplit;
import stelarx.dp.DPTable;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;
import stelarx.tree.TreeNode;
import stelarx.tree.TreeParser;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/** Checks duplicate-invariant S1 cluster and split deduplication. */
public final class DuplicateAwareCandidateTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("expected work directory");
        Path work = Path.of(args[0]).toAbsolutePath();
        Files.createDirectories(work);
        Path input = work.resolve("tagged-multicopy.tre");
        Files.writeString(input,
            "((((A,A)D,B),(C,D)),(E,F));\n(((A,B),(C,D)),(E,F));\n",
            StandardCharsets.UTF_8);

        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry, false);
        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 31L);
        PrefixHashArrays pref = new PrefixHashArrays(trees, hasher);

        ClusterTable oldClusters = new ClusterTable(trees, pref, registry.size());
        DPTable oldDP = new DPTable(trees, pref, oldClusters);

        UniqueTaxonSubtreeHashes unique = new UniqueTaxonSubtreeHashes(trees, hasher);
        for (Tree tree : trees) {
            verifySpeciationHashes(tree.root, tree, unique, hasher);
            verifyOnlySpeciationNodesIndexed(tree.root, tree.treeIndex, unique);
        }
        ClusterTable newClusters = new ClusterTable(trees, pref, registry.size(), unique);
        DPTable newDP = new DPTable(trees, pref, newClusters, unique);
        PartitionTable newPartitions = new PartitionTable(trees, pref, unique);

        check(oldClusters.size() == 12, "occurrence-sensitive cluster baseline");
        check(newClusters.size() == 10, "duplicate-invariant cluster count");
        check(oldDP.numUniqueSplits() == 8, "occurrence-sensitive split baseline");
        check(newDP.numUniqueSplits() == 5, "duplicate-invariant split count");
        check(newDP.numOverlappingSpeciationNodesSkipped() == 0,
            "correctly tagged speciation was rejected");
        check(newPartitions.entries().stream()
                .allMatch(entry -> entry.exemplar.size1 == entry.exemplar.hash1.size
                    && entry.exemplar.size2 == entry.exemplar.hash2.size
                    && entry.exemplar.size3 == entry.exemplar.hash3.size),
            "partition exemplar retained occurrence-sensitive sizes");
        check(newDP.getRootHash().size == registry.size(), "root distinct-taxon size");

        Tree first = trees.get(0);
        TreeNode duplication = first.root.left.left.left;
        check(!unique.contains(0, duplication),
            "duplication-rooted subtree received a standalone hash");
        BitSet singletonMembers = new BitSet(registry.size());
        singletonMembers.set(duplication.left.taxonId);
        ClusterHash singletonA = hashMembers(singletonMembers, hasher);
        check(newClusters.get(singletonA).frequency == 4,
            "duplication child set was not retained exactly once as a valid split side");

        TreeNode aabSpeciation = first.root.left.left;
        ClusterHash ab = unique.get(0, aabSpeciation);
        check(ab.size == 2, "speciation subtree did not deduplicate A");
        check(newClusters.get(ab).exemplar.size == 2,
            "cluster exemplar retained occurrence-sensitive size");

        for (var entry : newDP.entries()) {
            for (BipartitionSplit split : entry.getValue()) {
                check(split.lo.size + split.hi.size == entry.getKey().size,
                    "DP split contains overlapping species");
                check(newClusters.get(split.lo) != null
                        && newClusters.get(split.hi) != null,
                    "speciation transition side has no cluster exemplar");
            }
        }

        verifyOverlappingComplement(work);

        System.out.println("STELAR-Pro duplicate-aware S1 candidates: PASS");
    }

    /** A species present inside and outside sub(u) must remain in its complement. */
    private static void verifyOverlappingComplement(Path work) throws Exception {
        Path input = work.resolve("overlapping-complement.tre");
        Files.writeString(input, "(((A,B),C),(A,D))D;\n", StandardCharsets.UTF_8);
        TaxonRegistry registry = new TaxonRegistry();
        List<Tree> trees = TreeParser.parseGeneTrees(input.toString(), registry, false);
        TaxonHasher hasher = new TaxonHasher(registry.size(), 2, 41L);
        UniqueTaxonSubtreeHashes unique = new UniqueTaxonSubtreeHashes(trees, hasher);
        Tree tree = trees.get(0);

        verifySpeciationHashes(tree.root, tree, unique, hasher);
        verifyOnlySpeciationNodesIndexed(tree.root, tree.treeIndex, unique);

        TreeNode ab = tree.root.left.left;
        ClusterHash outsideAb = unique.getComplement(0, ab);
        check(outsideAb.size == 3,
            "outside hash removed a species that also has an outside copy");
    }

    /** Independent set oracle for each speciation-rooted bipartition record. */
    private static BitSet verifySpeciationHashes(
            TreeNode node, Tree tree, UniqueTaxonSubtreeHashes unique,
            TaxonHasher hasher) {
        BitSet members = new BitSet(hasher.numTaxa());
        if (node.isLeaf()) {
            members.set(node.taxonId);
        } else if (node.isPolytomous()) {
            for (int childIndex = 0; childIndex < node.children.length; childIndex++) {
                BitSet childMembers = verifySpeciationHashes(
                    node.children[childIndex], tree, unique, hasher);
                if (node.isSpeciation()) {
                    verifyHash(unique.getChild(tree.treeIndex, node, childIndex),
                        childMembers, hasher, "speciation child");
                }
                members.or(childMembers);
            }
        } else {
            BitSet left = verifySpeciationHashes(node.left, tree, unique, hasher);
            BitSet right = verifySpeciationHashes(node.right, tree, unique, hasher);
            if (node.isSpeciation()) {
                verifyHash(unique.getChild(tree.treeIndex, node, 0), left,
                    hasher, "left speciation child");
                verifyHash(unique.getChild(tree.treeIndex, node, 1), right,
                    hasher, "right speciation child");
            }
            members.or(left);
            members.or(right);
        }

        if (node.isSpeciation()) {
            verifyHash(unique.get(tree.treeIndex, node), members, hasher,
                "speciation subtree");
            BitSet outside = new BitSet(hasher.numTaxa());
            for (int position = 0; position < tree.leafCount; position++) {
                if (position < node.rangeStart || position >= node.rangeEnd) {
                    outside.set(tree.postorderArray[position]);
                }
            }
            verifyHash(unique.getComplement(tree.treeIndex, node), outside, hasher,
                "outside-subtree");
        }
        return members;
    }

    private static void verifyOnlySpeciationNodesIndexed(
            TreeNode node, int treeIndex, UniqueTaxonSubtreeHashes unique) {
        check(unique.contains(treeIndex, node) == node.isSpeciation(),
            "hash index contains a non-speciation node");
        if (node.isLeaf()) return;
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) {
                verifyOnlySpeciationNodesIndexed(child, treeIndex, unique);
            }
        } else {
            verifyOnlySpeciationNodesIndexed(node.left, treeIndex, unique);
            verifyOnlySpeciationNodesIndexed(node.right, treeIndex, unique);
        }
    }

    private static void verifyHash(ClusterHash actual, BitSet members,
                                   TaxonHasher hasher, String description) {
        ClusterHash expected = hashMembers(members, hasher);
        check(actual.size == expected.size, description + " size mismatch");
        check(Arrays.equals(actual.sums, expected.sums),
            description + " sum hash mismatch");
        check(Arrays.equals(actual.xors, expected.xors),
            description + " XOR hash mismatch");
    }

    private static ClusterHash hashMembers(BitSet members, TaxonHasher hasher) {
        long[] sums = new long[hasher.numSeeds()];
        long[] xors = new long[hasher.numSeeds()];
        for (int taxon = members.nextSetBit(0); taxon >= 0;
                taxon = members.nextSetBit(taxon + 1)) {
            for (int seed = 0; seed < hasher.numSeeds(); seed++) {
                long value = hasher.get(seed, taxon);
                sums[seed] += value;
                xors[seed] ^= value;
            }
        }
        return new ClusterHash(sums, xors, members.cardinality(), hasher.numSeeds());
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
