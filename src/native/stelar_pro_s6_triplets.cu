// S6 CUDA support builder: independent outgroups occupy independent blocks;
// each thread streams distinct gene partitions and atomically accumulates the
// exact number of speciation-rooted ab|c observations.
#include <jni.h>
#include <cuda_runtime.h>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(cudaError_t result, const char* operation) {
    if (result != cudaSuccess)
        throw std::runtime_error(std::string("S6 CUDA ") + operation + ": "
            + cudaGetErrorString(result));
}

__device__ bool contains(const int* taxa, int lo, int hi, int needle) {
    const int end = hi;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (taxa[mid] < needle) lo = mid + 1;
        else hi = mid;
    }
    return lo < end && taxa[lo] == needle;
}

__global__ void buildTensor(int n, int partCount, const int* parts,
                            const int* children, const int* taxa,
                            unsigned long long* tensor) {
    const int c = blockIdx.x;
    const size_t base = static_cast<size_t>(c) * n * n;
    for (int p = threadIdx.x; p < partCount; p += blockDim.x) {
        for (int out = parts[p]; out < parts[p + 1]; out++) {
            if (!contains(taxa, children[out], children[out + 1], c)) continue;
            for (int side = parts[p]; side < parts[p + 1]; side++) {
                if (side == out) continue;
                const int end = children[side + 1];
                for (int i = children[side]; i < end; i++) {
                    const int a = taxa[i];
                    if (a == c) continue;
                    for (int j = i + 1; j < end; j++) {
                        const int b = taxa[j];
                        if (b != c) atomicAdd(tensor + base + a * n + b, 1ULL);
                    }
                }
            }
        }
    }
}

}  // namespace

extern "C" JNIEXPORT jlongArray JNICALL
Java_stelarx_pro_S6TripletSupport_buildNativeGPU(JNIEnv* env, jclass, jint n,
    jintArray jPartOffsets, jintArray jChildOffsets, jintArray jTaxa) {
    int *dParts = nullptr, *dChildren = nullptr, *dTaxa = nullptr;
    unsigned long long* dTensor = nullptr;
    auto release = [&] {
        cudaFree(dParts); cudaFree(dChildren); cudaFree(dTaxa); cudaFree(dTensor);
    };
    try {
        const int partLength = env->GetArrayLength(jPartOffsets);
        const int childLength = env->GetArrayLength(jChildOffsets);
        const int taxonLength = env->GetArrayLength(jTaxa);
        const int64_t cells = int64_t(n) * n * n;
        if (n < 3 || partLength < 2 || cells > 2147483647)
            throw std::invalid_argument("invalid S6 CUDA tensor dimensions");
        std::vector<int> parts(partLength), children(childLength), taxa(taxonLength);
        env->GetIntArrayRegion(jPartOffsets, 0, partLength, parts.data());
        env->GetIntArrayRegion(jChildOffsets, 0, childLength, children.data());
        env->GetIntArrayRegion(jTaxa, 0, taxonLength, taxa.data());
        if (env->ExceptionCheck()) return nullptr;
        check(cudaMalloc(&dParts, partLength * sizeof(int)), "allocating parts");
        check(cudaMalloc(&dChildren, childLength * sizeof(int)), "allocating children");
        check(cudaMalloc(&dTaxa, taxonLength * sizeof(int)), "allocating taxa");
        check(cudaMalloc(&dTensor, static_cast<size_t>(cells) * sizeof(uint64_t)),
            "allocating tensor");
        check(cudaMemcpy(dParts, parts.data(), partLength * sizeof(int),
            cudaMemcpyHostToDevice), "copying parts");
        check(cudaMemcpy(dChildren, children.data(), childLength * sizeof(int),
            cudaMemcpyHostToDevice), "copying children");
        check(cudaMemcpy(dTaxa, taxa.data(), taxonLength * sizeof(int),
            cudaMemcpyHostToDevice), "copying taxa");
        check(cudaMemset(dTensor, 0, static_cast<size_t>(cells) * sizeof(uint64_t)),
            "clearing tensor");
        buildTensor<<<n, 256>>>(n, partLength - 1, dParts, dChildren, dTaxa,
            dTensor);
        check(cudaGetLastError(), "launching support kernel");
        check(cudaDeviceSynchronize(), "synchronizing support kernel");
        std::vector<jlong> result(static_cast<size_t>(cells));
        check(cudaMemcpy(result.data(), dTensor,
            static_cast<size_t>(cells) * sizeof(uint64_t), cudaMemcpyDeviceToHost),
            "copying support tensor");
        release();
        jlongArray output = env->NewLongArray(static_cast<jsize>(cells));
        if (output != nullptr)
            env->SetLongArrayRegion(output, 0, static_cast<jsize>(cells),
                result.data());
        return output;
    } catch (const std::exception& error) {
        release();
        jclass cls = env->FindClass("java/lang/IllegalStateException");
        if (cls) env->ThrowNew(cls, error.what());
        return nullptr;
    }
}
