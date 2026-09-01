package stelarx;

import stelarx.cluster.ClusterHash;
import stelarx.hash.PrefixHashArrays;
import stelarx.hash.TaxonHasher;
import stelarx.partition.Partition;
import stelarx.partition.PartitionTable;
import stelarx.taxon.TaxonRegistry;
import stelarx.tree.Tree;

import java.io.*;
import java.util.*;

/**
 * Verifies Phase-4 rooted child-partition extraction.
 *
 * Key checks:
 *   1. Every part's stored size/hash matches its distinct species set.
 *   2. Root partitions are retained (size3 == 0 is valid).
 *   3. Total frequency equals the number of eligible speciation nodes.
 *   4. For small inputs: show taxa in each part.
 */
public class Phase4Verifier {

    public static void dump(List<Tree> trees, TaxonRegistry registry,
                            TaxonHasher hasher, PrefixHashArrays pref,
                            PartitionTable partTable,
                            String outFile) throws IOException {
        PrintStream out = (outFile != null)
            ? new PrintStream(new FileOutputStream(outFile)) : System.out;

        int n = registry.size();
        int k = trees.size();
        int m = pref.numSeeds();

        out.printf("=== Phase 4 Rooted Child-Partition Verification ===%n");
        out.printf("Taxa: %d  Trees: %d  Seeds: %d%n", n, k, m);
        out.printf("Unique rooted child partitions: %d%n%n", partTable.size());

        // Expected eligible speciation nodes, root included for binary trees.
        int expectedTotal = 0;
        for (Tree t : trees) expectedTotal += countEligiblePartitions(t.root);
        out.printf("Expected speciation-rooted partitions: %d%n%n", expectedTotal);

        int fails = 0;

        // Check every part against a direct distinct-species oracle. Part sets
        // can overlap in multicopy trees, so their sizes need not sum to the
        // number of leaves or even to the number of distinct taxa in the tree.
        for (PartitionTable.Entry e : partTable.entries()) {
            Partition p = e.exemplar;
            Tree exemplarTree = trees.get(p.treeIndex);
            int childCount = p.d - 1;
            for (int part = 0; part < p.d; part++) {
                BitSet members;
                if (part < childCount) {
                    int start = p.d == 3 ? (part == 0 ? p.leftStart : p.rightStart)
                        : p.partStarts[part];
                    int end = p.d == 3 ? (part == 0 ? p.leftEnd : p.rightEnd)
                        : p.partEnds[part];
                    members = rangeMembers(exemplarTree, start, end, false, n);
                } else {
                    int start = p.d == 3 ? p.leftStart : p.partStarts[0];
                    int end = p.d == 3 ? p.rightEnd : p.partEnds[childCount - 1];
                    members = rangeMembers(exemplarTree, start, end, true, n);
                }
                ClusterHash actualHash = partitionHash(p, part);
                int actualSize = partitionSize(p, part);
                ClusterHash expectedHash = hashMembers(members, hasher);
                if (actualSize != members.cardinality()
                        || !actualHash.equals(expectedHash)) {
                    out.printf("FAIL: part %d distinct-set mismatch in %s%n", part, p);
                    fails++;
                }
            }
        }
        int observedTotal = partTable.entries().stream().mapToInt(e -> e.frequency).sum();
        if (observedTotal != expectedTotal) {
            out.printf("FAIL: partition frequency total=%d, expected=%d%n",
                observedTotal, expectedTotal);
            fails++;
        }

        out.printf("%n--- Summary ---%n");
        if (fails == 0) out.println("ALL ASSERTIONS PASSED");
        else            out.printf("%d FAILURES%n", fails);

        // Size distribution of part3 (tells us how "balanced" tripartitions are)
        out.printf("%n--- part3 size distribution (complement size) ---%n");
        Map<Integer, Integer> dist = new TreeMap<>();
        for (PartitionTable.Entry e : partTable.entries()) {
            dist.merge(e.exemplar.size3, 1, Integer::sum);
        }
        dist.forEach((sz, cnt) -> out.printf("  part3_size=%3d : %d tripartitions%n", sz, cnt));

        // Small-input: show all tripartitions with taxon names
        if (n <= 8) {
            out.printf("%n--- All tripartitions (small input) ---%n");
            List<PartitionTable.Entry> sorted = new ArrayList<>(partTable.entries());
            sorted.sort(Comparator.comparingInt(e -> e.exemplar.size1));
            for (PartitionTable.Entry e : sorted) {
                Partition p = e.exemplar;
                Tree t = trees.get(p.treeIndex);
                String s1 = rangeNames(t, p.leftStart,  p.leftEnd,  false, registry);
                String s2 = rangeNames(t, p.rightStart, p.rightEnd, false, registry);
                String s3 = rangeNames(t, p.leftStart,  p.rightEnd, true,  registry);
                out.printf("  freq=%d  {%s} | {%s} | {%s}%n", e.frequency, s1, s2, s3);
            }
        }

        if (outFile != null) out.close();
    }

    private static String rangeNames(Tree tree, int lo, int hi, boolean complement,
                                     TaxonRegistry registry) {
        StringBuilder sb = new StringBuilder();
        if (!complement) {
            for (int i = lo; i < hi; i++) {
                if (sb.length() > 0) sb.append(",");
                sb.append(registry.getName(tree.postorderArray[i]));
            }
        } else {
            for (int i = 0; i < tree.leafCount; i++) {
                if (i >= lo && i < hi) continue;
                if (sb.length() > 0) sb.append(",");
                sb.append(registry.getName(tree.postorderArray[i]));
            }
        }
        return sb.toString();
    }

    private static BitSet rangeMembers(Tree tree, int lo, int hi,
                                       boolean complement, int numTaxa) {
        BitSet members = new BitSet(numTaxa);
        for (int position = 0; position < tree.leafCount; position++) {
            if (complement == (position >= lo && position < hi)) continue;
            members.set(tree.postorderArray[position]);
        }
        return members;
    }

    private static ClusterHash hashMembers(BitSet members, TaxonHasher hasher) {
        int m = hasher.numSeeds();
        long[] sums = new long[m];
        long[] xors = new long[m];
        for (int taxon = members.nextSetBit(0); taxon >= 0;
                taxon = members.nextSetBit(taxon + 1)) {
            for (int seed = 0; seed < m; seed++) {
                long value = hasher.get(seed, taxon);
                sums[seed] += value;
                xors[seed] ^= value;
            }
        }
        return new ClusterHash(sums, xors, members.cardinality(), m);
    }

    private static ClusterHash partitionHash(Partition partition, int part) {
        if (partition.d != 3) return partition.hashes[part];
        return switch (part) {
            case 0 -> partition.hash1;
            case 1 -> partition.hash2;
            case 2 -> partition.hash3;
            default -> throw new IllegalArgumentException("invalid partition part");
        };
    }

    private static int partitionSize(Partition partition, int part) {
        if (partition.d != 3) return partition.sizes[part];
        return switch (part) {
            case 0 -> partition.size1;
            case 1 -> partition.size2;
            case 2 -> partition.size3;
            default -> throw new IllegalArgumentException("invalid partition part");
        };
    }

    private static int countEligiblePartitions(stelarx.tree.TreeNode node) {
        if (node.isLeaf()) return 0;
        int count = node.isSpeciation() && (!node.isPolytomous() || !node.isRoot()) ? 1 : 0;
        if (node.isPolytomous()) {
            for (var child : node.children) count += countEligiblePartitions(child);
        } else {
            count += countEligiblePartitions(node.left);
            count += countEligiblePartitions(node.right);
        }
        return count;
    }
}
