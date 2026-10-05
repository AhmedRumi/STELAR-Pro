#ifndef STELAR_PRO_S5_CUDA_API_H
#define STELAR_PRO_S5_CUDA_API_H

// Plain C interface keeps nvcc away from ASTRAL-Pro's large C++ search unit.
struct S5CudaNode {
    int up, small, large, dup;
    long long x, y, z, tx, ty, tz, q;
};

struct S5CudaOperation {
    int color;  // -2 requests the current quartet score; -1 removes a taxon
    int taxon;
    int scoreSlot;
};

extern "C" {
void* s5CudaCreate(int genes, int taxa, int nodeCount, const S5CudaNode* nodes,
                   const int* leafOffsets, int leafNodeCount,
                   const int* leafNodes);
void s5CudaExecute(void* context, int operationCount,
                   const S5CudaOperation* operations, int scoreCount,
                   long long* scoreSums);
void s5CudaDestroy(void* context);
}

#endif
