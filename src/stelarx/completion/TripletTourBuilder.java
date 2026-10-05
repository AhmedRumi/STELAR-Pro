package stelarx.completion;

import stelarx.tree.Tree;
import stelarx.tree.TreeNode;

import java.util.Arrays;
import java.util.IdentityHashMap;

/**
 * Euler/RMQ data needed by the rooted-triplet similarity matrix.
 *
 * <p>Unlike {@link EulerTourBuilder#buildFull}, this builder stores only the
 * descendant-leaf count of the node at each Euler position. It deliberately
 * does not construct the quartet bridge payloads (F and child-of-LCA data).
 */
final class TripletTourBuilder {
    private TripletTourBuilder() {}

    static final class CompactData {
        final short[] depths;
        final short[] nodeSizes;
        final char[][] sparseArgmin;
        final int[] firstOcc;
        final int tourLen;
        final int log;
        final int leafCount;

        CompactData(short[] depths, short[] nodeSizes, char[][] sparseArgmin,
                    int[] firstOcc, int tourLen, int log, int leafCount) {
            this.depths = depths;
            this.nodeSizes = nodeSizes;
            this.sparseArgmin = sparseArgmin;
            this.firstOcc = firstOcc;
            this.tourLen = tourLen;
            this.log = log;
            this.leafCount = leafCount;
        }
    }

    static final class WideData {
        final int[] depths;
        final int[] nodeSizes;
        final byte[][] microArgmin;
        final int[][] macroArgmin;
        final int[] firstOcc;
        final int tourLen;
        final int blockCount;
        final int macroLog;
        final int leafCount;

        WideData(int[] depths, int[] nodeSizes, byte[][] microArgmin,
                 int[][] macroArgmin, int[] firstOcc, int tourLen,
                 int blockCount, int macroLog, int leafCount) {
            this.depths = depths;
            this.nodeSizes = nodeSizes;
            this.microArgmin = microArgmin;
            this.macroArgmin = macroArgmin;
            this.firstOcc = firstOcc;
            this.tourLen = tourLen;
            this.blockCount = blockCount;
            this.macroLog = macroLog;
            this.leafCount = leafCount;
        }
    }

    static CompactData buildCompact(Tree tree, int n) {
        requireSingleCopy(tree);
        int len = EulerTourBuilder.tourLength(tree);
        if (len > Character.MAX_VALUE + 1) {
            throw new IllegalArgumentException("Triplet Euler tour has " + len
                + " positions; compact RMQ supports at most "
                + (Character.MAX_VALUE + 1));
        }

        short[] depths = new short[len];
        short[] nodeSizes = new short[len];
        int[] firstOcc = new int[n];
        Arrays.fill(firstOcc, -1);
        IdentityHashMap<TreeNode, Integer> sizes = subtreeSizes(tree.root);
        int[] cursor = {0};
        emitCompact(tree.root, 0, sizes, depths, nodeSizes, firstOcc, cursor);

        int log = sparseLog(len);
        char[][] argmin = new char[log][len];
        for (int p = 0; p < len; p++) argmin[0][p] = (char)p;
        for (int level = 1; level < log; level++) {
            int half = 1 << (level - 1);
            int end = len - (1 << level) + 1;
            for (int p = 0; p < end; p++) {
                int left = argmin[level - 1][p];
                int right = argmin[level - 1][p + half];
                argmin[level][p] = (char)(depths[left] <= depths[right] ? left : right);
            }
        }
        return new CompactData(depths, nodeSizes, argmin, firstOcc,
            len, log, tree.leafCount);
    }

    static WideData buildWide(Tree tree, int n) {
        requireSingleCopy(tree);
        int len = EulerTourBuilder.tourLength(tree);
        int[] depths = new int[len];
        int[] nodeSizes = new int[len];
        int[] firstOcc = new int[n];
        Arrays.fill(firstOcc, -1);
        IdentityHashMap<TreeNode, Integer> sizes = subtreeSizes(tree.root);
        int[] cursor = {0};
        emitWide(tree.root, 0, sizes, depths, nodeSizes, firstOcc, cursor);

        int microLog = EulerTourBuilder.WIDE_MICRO_LOG;
        int blockSize = EulerTourBuilder.WIDE_BLOCK_SIZE;
        byte[][] micro = new byte[microLog][len];
        for (int p = 0; p < len; p++) micro[0][p] = (byte)(p & 0xff);
        for (int level = 1; level < microLog; level++) {
            int width = 1 << level;
            int half = width >>> 1;
            for (int start = 0; start < len; start += blockSize) {
                int end = Math.min(start + blockSize, len);
                for (int p = start; p <= end - width; p++) {
                    int left = start + Byte.toUnsignedInt(micro[level - 1][p]);
                    int right = start + Byte.toUnsignedInt(micro[level - 1][p + half]);
                    micro[level][p] = (byte)((depths[left] <= depths[right]
                        ? left : right) - start);
                }
            }
        }

        int blocks = (len + blockSize - 1) / blockSize;
        int macroLog = sparseLog(blocks);
        int[][] macro = new int[macroLog][blocks];
        for (int block = 0; block < blocks; block++) {
            int lo = block * blockSize;
            int hi = Math.min(lo + blockSize, len);
            int best = lo;
            for (int p = lo + 1; p < hi; p++) {
                if (depths[p] < depths[best]) best = p;
            }
            macro[0][block] = best;
        }
        for (int level = 1; level < macroLog; level++) {
            int half = 1 << (level - 1);
            int end = blocks - (1 << level) + 1;
            for (int block = 0; block < end; block++) {
                int left = macro[level - 1][block];
                int right = macro[level - 1][block + half];
                macro[level][block] = depths[left] <= depths[right] ? left : right;
            }
        }
        return new WideData(depths, nodeSizes, micro, macro, firstOcc,
            len, blocks, macroLog, tree.leafCount);
    }

    static int queryCompact(CompactData data, int lo, int hi) {
        int width = hi - lo + 1;
        int level = 31 - Integer.numberOfLeadingZeros(width);
        int rightStart = hi - (1 << level) + 1;
        int left = data.sparseArgmin[level][lo];
        int right = data.sparseArgmin[level][rightStart];
        return data.depths[left] <= data.depths[right] ? left : right;
    }

    static int queryWide(WideData data, int lo, int hi) {
        int blockSize = EulerTourBuilder.WIDE_BLOCK_SIZE;
        int leftBlock = lo / blockSize;
        int rightBlock = hi / blockSize;
        if (leftBlock == rightBlock) return queryMicro(data, lo, hi);

        int best = queryMicro(data, lo, (leftBlock + 1) * blockSize - 1);
        int firstWhole = leftBlock + 1;
        int lastWhole = rightBlock - 1;
        if (firstWhole <= lastWhole) {
            int count = lastWhole - firstWhole + 1;
            int level = 31 - Integer.numberOfLeadingZeros(count);
            int second = lastWhole - (1 << level) + 1;
            int left = data.macroArgmin[level][firstWhole];
            int right = data.macroArgmin[level][second];
            int middle = data.depths[left] <= data.depths[right] ? left : right;
            if (data.depths[middle] < data.depths[best]) best = middle;
        }
        int right = queryMicro(data, rightBlock * blockSize, hi);
        return data.depths[right] < data.depths[best] ? right : best;
    }

    private static int queryMicro(WideData data, int lo, int hi) {
        int width = hi - lo + 1;
        int level = 31 - Integer.numberOfLeadingZeros(width);
        int second = hi - (1 << level) + 1;
        int start = (lo / EulerTourBuilder.WIDE_BLOCK_SIZE)
            * EulerTourBuilder.WIDE_BLOCK_SIZE;
        int left = start + Byte.toUnsignedInt(data.microArgmin[level][lo]);
        int right = start + Byte.toUnsignedInt(data.microArgmin[level][second]);
        return data.depths[left] <= data.depths[right] ? left : right;
    }

    private static IdentityHashMap<TreeNode, Integer> subtreeSizes(TreeNode root) {
        IdentityHashMap<TreeNode, Integer> sizes = new IdentityHashMap<>();
        computeSize(root, sizes);
        return sizes;
    }

    private static int computeSize(TreeNode node, IdentityHashMap<TreeNode, Integer> sizes) {
        int size;
        if (node.isLeaf()) {
            size = 1;
        } else if (node.isPolytomous()) {
            size = 0;
            for (TreeNode child : node.children) size += computeSize(child, sizes);
        } else {
            size = computeSize(node.left, sizes) + computeSize(node.right, sizes);
        }
        sizes.put(node, size);
        return size;
    }

    private static void emitCompact(TreeNode node, int depth,
                                    IdentityHashMap<TreeNode, Integer> sizes,
                                    short[] depths, short[] nodeSizes,
                                    int[] firstOcc, int[] cursor) {
        int pos = cursor[0]++;
        depths[pos] = (short)depth;
        nodeSizes[pos] = (short)(int)sizes.get(node);
        if (node.isLeaf()) {
            firstOcc[node.taxonId] = pos;
            return;
        }
        TreeNode[] children = children(node);
        emitCompact(children[0], depth + 1, sizes, depths, nodeSizes, firstOcc, cursor);
        for (int i = 1; i < children.length; i++) {
            int middle = cursor[0]++;
            depths[middle] = (short)depth;
            nodeSizes[middle] = (short)(int)sizes.get(node);
            emitCompact(children[i], depth + 1, sizes, depths, nodeSizes, firstOcc, cursor);
        }
    }

    private static void emitWide(TreeNode node, int depth,
                                 IdentityHashMap<TreeNode, Integer> sizes,
                                 int[] depths, int[] nodeSizes,
                                 int[] firstOcc, int[] cursor) {
        int pos = cursor[0]++;
        depths[pos] = depth;
        nodeSizes[pos] = sizes.get(node);
        if (node.isLeaf()) {
            firstOcc[node.taxonId] = pos;
            return;
        }
        TreeNode[] children = children(node);
        emitWide(children[0], depth + 1, sizes, depths, nodeSizes, firstOcc, cursor);
        for (int i = 1; i < children.length; i++) {
            int middle = cursor[0]++;
            depths[middle] = depth;
            nodeSizes[middle] = sizes.get(node);
            emitWide(children[i], depth + 1, sizes, depths, nodeSizes, firstOcc, cursor);
        }
    }

    private static TreeNode[] children(TreeNode node) {
        return node.isPolytomous()
            ? node.children
            : new TreeNode[] { node.left, node.right };
    }

    private static int sparseLog(int value) {
        int log = 1;
        while ((1 << log) <= value && log < 31) log++;
        return log;
    }

    private static void requireSingleCopy(Tree tree) {
        if (tree.distinctTaxonCount != tree.leafCount) {
            throw new IllegalArgumentException(
                "Triplet similarity requires single-copy DISCO trees");
        }
    }
}
