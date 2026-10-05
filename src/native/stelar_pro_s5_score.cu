// Per-locus A-Pro quartet counters. A CUDA worker owns one gene family and
// applies the same ordered color updates as multitree.hpp::Partition::update.
#include "stelar_pro_s5_cuda_api.h"

#include <cuda_runtime.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(cudaError_t status, const char* action) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("S5 CUDA ") + action + ": "
                                 + cudaGetErrorString(status));
    }
}

__device__ inline long long& countForColor(S5CudaNode& n, int color) {
    return color == 0 ? n.z : color == 1 ? n.x : n.y;
}

__device__ inline void updateDuplication(S5CudaNode& n, S5CudaNode* nodes) {
    const S5CudaNode& u = nodes[n.small];
    const S5CudaNode& v = nodes[n.large];
    n.tx = u.tx + v.tx;
    n.ty = u.ty + v.ty;
    n.tz = u.tz + v.tz;
}

__device__ inline long long updateSpeciation(S5CudaNode& n,
                                              S5CudaNode* nodes) {
    const S5CudaNode& u = nodes[n.small];
    const S5CudaNode& v = nodes[n.large];
    n.tx = u.tx + v.tx + u.y * v.z * (u.y + v.z - 2)
                       + u.z * v.y * (u.z + v.y - 2);
    n.ty = u.ty + v.ty + u.x * v.z * (u.x + v.z - 2)
                       + u.z * v.x * (u.z + v.x - 2);
    n.tz = u.tz + v.tz + u.x * v.y * (u.x + v.y - 2)
                       + u.y * v.x * (u.y + v.x - 2);
    const long long old = n.q;
    n.q = u.x * v.tx + u.y * v.ty + u.z * v.tz
        + v.x * u.tx + v.y * u.ty + v.z * u.tz
        + (u.x - 1) * u.x * v.y * v.z
        + (u.y - 1) * u.y * v.x * v.z
        + (u.z - 1) * u.z * v.x * v.y;
    return n.q - old;
}

__global__ void applyOperations(S5CudaNode* nodes, const int* offsets,
                                 const int* leafNodes, int* colors,
                                 long long* totals,
                                 const S5CudaOperation* operations,
                                 int operationCount, int genes, int taxa,
                                 long long* perGeneScores) {
    const int gene = blockIdx.x * blockDim.x + threadIdx.x;
    if (gene >= genes) return;
    long long total = totals[gene];
    for (int k = 0; k < operationCount; ++k) {
        const S5CudaOperation op = operations[k];
        if (op.color == -2) {
            perGeneScores[(size_t)op.scoreSlot * genes + gene] = total / 2;
            continue;
        }
        const int row = gene * taxa + op.taxon;
        const int old = colors[row];
        if (old == op.color) continue;
        for (int j = offsets[row]; j < offsets[row + 1]; ++j) {
            int u = leafNodes[j];
            if (old != -1) --countForColor(nodes[u], old);
            if (op.color != -1) ++countForColor(nodes[u], op.color);
            int w = nodes[u].up;
            while (w != -1 && (!nodes[w].dup || nodes[w].large == u)) {
                if (old != -1) --countForColor(nodes[w], old);
                if (op.color != -1) ++countForColor(nodes[w], op.color);
                if (nodes[w].dup) updateDuplication(nodes[w], nodes);
                else total += updateSpeciation(nodes[w], nodes);
                u = w;
                w = nodes[u].up;
            }
            if (w != -1) updateDuplication(nodes[w], nodes);
        }
        colors[row] = op.color;
    }
    totals[gene] = total;
}

__global__ void reduceScores(const long long* perGene, long long* sums,
                              int genes) {
    const int slot = blockIdx.x;
    __shared__ long long partial[256];
    long long sum = 0;
    for (int gene = threadIdx.x; gene < genes; gene += blockDim.x) {
        sum += perGene[(size_t)slot * genes + gene];
    }
    partial[threadIdx.x] = sum;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) partial[threadIdx.x] += partial[threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0) sums[slot] = partial[0];
}

struct Context {
    int genes, taxa;
    S5CudaNode* nodes = nullptr;
    int *offsets = nullptr, *leafNodes = nullptr, *colors = nullptr;
    long long* totals = nullptr;

    ~Context() {
        cudaFree(nodes); cudaFree(offsets); cudaFree(leafNodes);
        cudaFree(colors); cudaFree(totals);
    }
};

}  // namespace

extern "C" void* s5CudaCreate(int genes, int taxa, int nodeCount,
                                const S5CudaNode* nodes, const int* leafOffsets,
                                int leafNodeCount, const int* leafNodes) {
    if (genes < 1 || taxa < 4 || nodeCount < 1 || leafNodeCount < 1) {
        throw std::invalid_argument("Invalid S5 CUDA gene-tree dimensions");
    }
    int devices = 0;
    check(cudaGetDeviceCount(&devices), "querying devices");
    if (devices < 1) throw std::runtime_error("S5 GPU requested without a CUDA device");
    Context* ctx = new Context{genes, taxa};
    try {
        const size_t colorCount = static_cast<size_t>(genes) * taxa;
        check(cudaMalloc(&ctx->nodes, nodeCount * sizeof(S5CudaNode)), "allocating nodes");
        check(cudaMalloc(&ctx->offsets, (colorCount + 1) * sizeof(int)), "allocating offsets");
        check(cudaMalloc(&ctx->leafNodes, leafNodeCount * sizeof(int)), "allocating leaves");
        check(cudaMalloc(&ctx->colors, colorCount * sizeof(int)), "allocating colors");
        check(cudaMalloc(&ctx->totals, genes * sizeof(long long)), "allocating totals");
        check(cudaMemcpy(ctx->nodes, nodes, nodeCount * sizeof(S5CudaNode),
                         cudaMemcpyHostToDevice), "copying nodes");
        check(cudaMemcpy(ctx->offsets, leafOffsets, (colorCount + 1) * sizeof(int),
                         cudaMemcpyHostToDevice), "copying offsets");
        check(cudaMemcpy(ctx->leafNodes, leafNodes, leafNodeCount * sizeof(int),
                         cudaMemcpyHostToDevice), "copying leaves");
        check(cudaMemset(ctx->colors, 0xff, colorCount * sizeof(int)),
              "initializing colors");
        check(cudaMemset(ctx->totals, 0, genes * sizeof(long long)),
              "initializing totals");
        return ctx;
    } catch (...) {
        delete ctx;
        throw;
    }
}

extern "C" void s5CudaExecute(void* context, int operationCount,
                                const S5CudaOperation* operations, int scoreCount,
                                long long* scoreSums) {
    if (!context || operationCount < 1 || scoreCount < 0) {
        throw std::invalid_argument("Invalid S5 CUDA operation batch");
    }
    Context& ctx = *static_cast<Context*>(context);
    S5CudaOperation* dOps = nullptr;
    long long *dPerGene = nullptr, *dSums = nullptr;
    try {
        const size_t scoreCells = static_cast<size_t>(std::max(1, scoreCount)) * ctx.genes;
        check(cudaMalloc(&dOps, operationCount * sizeof(S5CudaOperation)),
              "allocating operations");
        check(cudaMalloc(&dPerGene, scoreCells * sizeof(long long)),
              "allocating per-gene scores");
        if (scoreCount > 0) {
            check(cudaMalloc(&dSums, scoreCount * sizeof(long long)),
                  "allocating score sums");
        }
        check(cudaMemcpy(dOps, operations, operationCount * sizeof(S5CudaOperation),
                         cudaMemcpyHostToDevice), "copying operations");
        applyOperations<<<(ctx.genes + 31) / 32, 32>>>(
            ctx.nodes, ctx.offsets, ctx.leafNodes, ctx.colors, ctx.totals,
            dOps, operationCount, ctx.genes, ctx.taxa, dPerGene);
        check(cudaGetLastError(), "launching color updates");
        if (scoreCount > 0) {
            reduceScores<<<scoreCount, 256>>>(dPerGene, dSums, ctx.genes);
            check(cudaGetLastError(), "launching score reduction");
            check(cudaMemcpy(scoreSums, dSums, scoreCount * sizeof(long long),
                             cudaMemcpyDeviceToHost), "copying quartet scores");
        } else {
            check(cudaDeviceSynchronize(), "finishing color updates");
        }
    } catch (...) {
        cudaFree(dOps); cudaFree(dPerGene); cudaFree(dSums);
        throw;
    }
    cudaFree(dOps); cudaFree(dPerGene); cudaFree(dSums);
}

extern "C" void s5CudaDestroy(void* context) {
    delete static_cast<Context*>(context);
}
