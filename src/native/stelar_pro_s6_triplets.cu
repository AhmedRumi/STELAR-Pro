// S6 CUDA support. The dense tensor builder remains for parity tests; the
// production path builds only the focused slices needed for insertions.
#include <jni.h>
#include <cuda_runtime.h>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <atomic>

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

// Insertion uses only triplets containing x. The inverted index builds that
// O(n²) slice by visiting child sides containing x, never the full tensor.
__global__ void focusedByMembership(int n, const int* parts,
                                    const int* children, const int* taxa,
                                    const int* partOfChild,
                                    const int* taxonChildren, int first, int last,
                                    int x, unsigned long long* output) {
    int membership = first + blockIdx.x * blockDim.x + threadIdx.x;
    if (membership >= last) return;
    int side = taxonChildren[membership];
    int p = partOfChild[side];
    unsigned long long* paired = output + static_cast<size_t>(n) * n;
    for (int out = parts[p]; out < parts[p + 1]; out++) {
        if (out == side) continue;
        int end = children[out + 1];
        for (int i = children[out]; i < end; i++) {
            int a = taxa[i];
            if (a == x) continue;
            for (int j = i + 1; j < end; j++) {
                int b = taxa[j];
                if (b != x) atomicAdd(output + a * n + b, 1ULL);
            }
        }
        for (int i = children[out]; i < end; i++) {
            int c = taxa[i];
            if (c == x) continue;
            for (int j = children[side]; j < children[side + 1]; j++) {
                int a = taxa[j];
                if (a != x && a != c) atomicAdd(paired + c * n + a, 1ULL);
            }
        }
    }
}

// The two-step search only asks about pairs in its fixed backbone. Filter
// child-side taxa before pair enumeration and emit a compact 2*m*m slice.
__global__ void focusedSubsetByMembership(int m, const int* parts,
                                          const int* children, const int* taxa,
                                          const int* partOfChild,
                                          const int* taxonChildren,
                                          const int* selectedMap,
                                          int first, int last,
                                          unsigned long long* output) {
    int membership = first + blockIdx.x * blockDim.x + threadIdx.x;
    if (membership >= last) return;
    int side = taxonChildren[membership];
    int p = partOfChild[side];
    unsigned long long* paired = output + static_cast<size_t>(m) * m;
    for (int out = parts[p]; out < parts[p + 1]; out++) {
        if (out == side) continue;
        for (int i = children[out]; i < children[out + 1]; i++) {
            int a = selectedMap[taxa[i]];
            if (a < 0) continue;
            for (int j = i + 1; j < children[out + 1]; j++) {
                int b = selectedMap[taxa[j]];
                if (b >= 0) atomicAdd(output + a * m + b, 1ULL);
            }
        }
        for (int i = children[out]; i < children[out + 1]; i++) {
            int c = selectedMap[taxa[i]];
            if (c < 0) continue;
            for (int j = children[side]; j < children[side + 1]; j++) {
                int a = selectedMap[taxa[j]];
                if (a >= 0 && a != c)
                    atomicAdd(paired + c * m + a, 1ULL);
            }
        }
    }
}

// For large gene sets, upload only the sides needed by this focused taxon.
// Taxa are already mapped to the selected backbone, so device memory is
// proportional to this request rather than all multicopy child memberships.
__global__ void focusedCompactByMembership(int m, const int* parts,
                                          const int* children, const int* taxa,
                                          const int* partOfChild,
                                          const int* memberships, int count,
                                          unsigned long long* output) {
    int membership = blockIdx.x * blockDim.x + threadIdx.x;
    if (membership >= count) return;
    int side = memberships[membership];
    int p = partOfChild[side];
    unsigned long long* paired = output + static_cast<size_t>(m) * m;
    for (int out = parts[p]; out < parts[p + 1]; out++) {
        if (out == side) continue;
        for (int i = children[out]; i < children[out + 1]; i++) {
            int a = taxa[i];
            for (int j = i + 1; j < children[out + 1]; j++)
                atomicAdd(output + a * m + taxa[j], 1ULL);
        }
        for (int i = children[out]; i < children[out + 1]; i++) {
            int c = taxa[i];
            for (int j = children[side]; j < children[side + 1]; j++) {
                int a = taxa[j];
                if (a != c) atomicAdd(paired + c * m + a, 1ULL);
            }
        }
    }
}

struct FocusIndex {
    int n, partCount;
    int *parts = nullptr, *children = nullptr, *taxa = nullptr;
    int *partOfChild = nullptr, *taxonChildren = nullptr;
    std::vector<int> offsets;
    // S5 does not keep a full guide index resident on the GPU. Match that
    // bounded-VRAM behavior for large S6 inputs while retaining CUDA scoring.
    bool staged = false;
    std::vector<int> hostParts, hostChildren, hostTaxa;
    std::vector<int> hostPartOfChild, hostTaxonChildren;
    std::vector<int> selected, selectedChildren, selectedTaxa;
    std::atomic_flag selectedBusy = ATOMIC_FLAG_INIT;
    ~FocusIndex() {
        cudaFree(parts); cudaFree(children); cudaFree(taxa);
        cudaFree(partOfChild); cudaFree(taxonChildren);
    }
};

struct CompactFocus {
    std::vector<int> parts{0}, children{0}, taxa, partOfChild, memberships;
};

struct SelectedGuard {
    std::atomic_flag& busy;
    explicit SelectedGuard(std::atomic_flag& flag) : busy(flag) {
        while (busy.test_and_set(std::memory_order_acquire)) {}
    }
    ~SelectedGuard() { busy.clear(std::memory_order_release); }
};

CompactFocus compactFocus(const FocusIndex& index, int x,
                          const std::vector<int>& childOffsets,
                          const std::vector<int>& childTaxa,
                          const std::vector<int>& map, bool mapped) {
    CompactFocus compact;
    // Java constructs each taxon's membership list in ascending child order,
    // hence partitions are contiguous here. Reuse a copied partition when x
    // appears in both children, without a huge per-focus remapping array.
    int lastPart = -1;
    int local = -1;
    for (int at = index.offsets[x]; at < index.offsets[x + 1]; at++) {
        int side = index.hostTaxonChildren[at];
        int p = index.hostPartOfChild[side];
        if (p != lastPart) {
            local = static_cast<int>(compact.parts.size()) - 1;
            lastPart = p;
            int firstSide = index.hostParts[p];
            int lastSide = index.hostParts[p + 1];
            for (int oldSide = firstSide; oldSide < lastSide; oldSide++) {
                for (int i = childOffsets[oldSide];
                     i < childOffsets[oldSide + 1]; i++) {
                    int selected = mapped ? childTaxa[i] : map[childTaxa[i]];
                    if (selected >= 0) compact.taxa.push_back(selected);
                }
                compact.children.push_back(static_cast<int>(compact.taxa.size()));
                compact.partOfChild.push_back(local);
            }
            compact.parts.push_back(static_cast<int>(compact.partOfChild.size()));
        }
        compact.memberships.push_back(compact.parts[local]
            + side - index.hostParts[p]);
    }
    return compact;
}

void prepareSelected(FocusIndex& index, const std::vector<jint>& selected,
                     const std::vector<int>& map) {
    if (index.selected.size() == selected.size()
            && std::equal(index.selected.begin(), index.selected.end(),
                selected.begin())) return;
    index.selected.assign(selected.begin(), selected.end());
    index.selectedChildren.clear();
    index.selectedTaxa.clear();
    index.selectedChildren.reserve(index.hostChildren.size());
    index.selectedChildren.push_back(0);
    for (size_t side = 0; side + 1 < index.hostChildren.size(); side++) {
        for (int at = index.hostChildren[side];
             at < index.hostChildren[side + 1]; at++) {
            int id = map[index.hostTaxa[at]];
            if (id >= 0) index.selectedTaxa.push_back(id);
        }
        index.selectedChildren.push_back(static_cast<int>(index.selectedTaxa.size()));
    }
}

void runCompactFocus(const FocusIndex& index, int x,
                     const std::vector<int>& childOffsets,
                     const std::vector<int>& childTaxa,
                     const std::vector<int>& map, bool mapped, int m,
                     unsigned long long* host) {
    CompactFocus compact = compactFocus(index, x, childOffsets,
        childTaxa, map, mapped);
    int *dParts = nullptr, *dChildren = nullptr, *dTaxa = nullptr;
    int *dPartOfChild = nullptr, *dMemberships = nullptr;
    unsigned long long* dOutput = nullptr;
    auto release = [&] {
        cudaFree(dParts); cudaFree(dChildren); cudaFree(dTaxa);
        cudaFree(dPartOfChild); cudaFree(dMemberships); cudaFree(dOutput);
    };
    try {
        auto upload = [](int** device, const std::vector<int>& source,
                         const char* label) {
            size_t bytes = std::max<size_t>(1, source.size()) * sizeof(int);
            check(cudaMalloc(device, bytes), label);
            if (!source.empty())
                check(cudaMemcpy(*device, source.data(), source.size() * sizeof(int),
                    cudaMemcpyHostToDevice), label);
        };
        upload(&dParts, compact.parts, "uploading focused parts");
        upload(&dChildren, compact.children, "uploading focused children");
        upload(&dTaxa, compact.taxa, "uploading focused taxa");
        upload(&dPartOfChild, compact.partOfChild, "uploading focused part lookup");
        upload(&dMemberships, compact.memberships, "uploading focused membership");
        size_t bytes = 2ULL * m * m * sizeof(uint64_t);
        check(cudaMalloc(&dOutput, bytes), "allocating compact focus output");
        check(cudaMemset(dOutput, 0, bytes), "clearing compact focus output");
        if (!compact.memberships.empty()) {
            int blocks = (compact.memberships.size() + 255) / 256;
            focusedCompactByMembership<<<blocks, 256>>>(m, dParts, dChildren,
                dTaxa, dPartOfChild, dMemberships,
                static_cast<int>(compact.memberships.size()), dOutput);
            check(cudaGetLastError(), "launching compact focus kernel");
        }
        check(cudaMemcpy(host, dOutput, bytes, cudaMemcpyDeviceToHost),
            "copying compact focus");
        release();
    } catch (...) {
        release();
        throw;
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

extern "C" JNIEXPORT jlong JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_createFocusIndex(JNIEnv* env, jclass,
    jint n, jintArray jParts, jintArray jChildren, jintArray jTaxa,
    jintArray jPartOfChild, jintArray jOffsets, jintArray jTaxonChildren) {
    FocusIndex* index = new FocusIndex();
    try {
        int pLen = env->GetArrayLength(jParts);
        int cLen = env->GetArrayLength(jChildren);
        int tLen = env->GetArrayLength(jTaxa);
        int mLen = env->GetArrayLength(jTaxonChildren);
        if (n < 3 || pLen < 2 || cLen < 3 ||
                env->GetArrayLength(jPartOfChild) != cLen - 1 ||
                env->GetArrayLength(jOffsets) != n + 1 || mLen != tLen)
            throw std::invalid_argument("invalid S6 focused index dimensions");
        index->n = n;
        index->partCount = pLen - 1;
        std::vector<int> p(pLen), c(cLen), t(tLen), member(mLen), partOf(cLen - 1);
        index->offsets.resize(n + 1);
        env->GetIntArrayRegion(jParts, 0, pLen, p.data());
        env->GetIntArrayRegion(jChildren, 0, cLen, c.data());
        env->GetIntArrayRegion(jTaxa, 0, tLen, t.data());
        env->GetIntArrayRegion(jPartOfChild, 0, cLen - 1, partOf.data());
        env->GetIntArrayRegion(jOffsets, 0, n + 1, index->offsets.data());
        env->GetIntArrayRegion(jTaxonChildren, 0, mLen, member.data());
        if (env->ExceptionCheck()) { delete index; return 0; }
        if (n > 256) {
            // Keep the large immutable index in host RAM. Each CUDA request
            // uploads only partitions containing its focused taxon and only
            // taxa from the current backbone subset.
            index->staged = true;
            index->hostParts = std::move(p);
            index->hostChildren = std::move(c);
            index->hostTaxa = std::move(t);
            index->hostPartOfChild = std::move(partOf);
            index->hostTaxonChildren = std::move(member);
            return reinterpret_cast<jlong>(index);
        }
        check(cudaMalloc(&index->parts, pLen * sizeof(int)), "allocating focused parts");
        check(cudaMalloc(&index->children, cLen * sizeof(int)), "allocating focused children");
        check(cudaMalloc(&index->taxa, tLen * sizeof(int)), "allocating focused taxa");
        check(cudaMalloc(&index->partOfChild, (cLen - 1) * sizeof(int)),
            "allocating focused part lookup");
        check(cudaMalloc(&index->taxonChildren, mLen * sizeof(int)),
            "allocating focused taxon lookup");
        check(cudaMemcpy(index->parts, p.data(), pLen * sizeof(int), cudaMemcpyHostToDevice),
            "copying focused parts");
        check(cudaMemcpy(index->children, c.data(), cLen * sizeof(int), cudaMemcpyHostToDevice),
            "copying focused children");
        check(cudaMemcpy(index->taxa, t.data(), tLen * sizeof(int), cudaMemcpyHostToDevice),
            "copying focused taxa");
        check(cudaMemcpy(index->partOfChild, partOf.data(), (cLen - 1) * sizeof(int),
            cudaMemcpyHostToDevice), "copying focused part lookup");
        check(cudaMemcpy(index->taxonChildren, member.data(), mLen * sizeof(int),
            cudaMemcpyHostToDevice), "copying focused taxon lookup");
        return reinterpret_cast<jlong>(index);
    } catch (const std::exception& error) {
        delete index;
        jclass cls = env->FindClass("java/lang/IllegalStateException");
        if (cls) env->ThrowNew(cls, error.what());
        return 0;
    }
}

extern "C" JNIEXPORT void JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_fillFocus(JNIEnv* env, jclass,
    jlong handle, jint x, jlongArray destination) {
    auto* index = reinterpret_cast<FocusIndex*>(handle);
    unsigned long long* output = nullptr;
    jlong* host = nullptr;
    try {
        if (!index || x < 0 || x >= index->n ||
                env->GetArrayLength(destination) != 2LL * index->n * index->n)
            throw std::invalid_argument("invalid S6 focused support request");
        if (index->staged) {
            std::vector<int> map(index->n);
            for (int taxon = 0; taxon < index->n; taxon++) map[taxon] = taxon;
            map[x] = -1;
            std::vector<jlong> result(2ULL * index->n * index->n);
            runCompactFocus(*index, x, index->hostChildren,
                index->hostTaxa, map, false, index->n,
                reinterpret_cast<unsigned long long*>(result.data()));
            env->SetLongArrayRegion(destination, 0,
                static_cast<jsize>(result.size()), result.data());
            return;
        }
        size_t bytes = 2ULL * index->n * index->n * sizeof(uint64_t);
        check(cudaMalloc(&output, bytes), "allocating focused output");
        check(cudaMemset(output, 0, bytes), "clearing focused output");
        int first = index->offsets[x], last = index->offsets[x + 1];
        if (last > first) {
            int blocks = (last - first + 255) / 256;
            focusedByMembership<<<blocks, 256>>>(index->n,
                index->parts, index->children, index->taxa,
                index->partOfChild, index->taxonChildren,
                first, last, x, output);
            check(cudaGetLastError(), "launching focused membership kernel");
        }
        host = env->GetLongArrayElements(destination, nullptr);
        if (!host) { cudaFree(output); return; }
        check(cudaMemcpy(host, output, bytes, cudaMemcpyDeviceToHost),
            "copying focused support");
        env->ReleaseLongArrayElements(destination, host, 0);
        cudaFree(output);
    } catch (const std::exception& error) {
        if (host) env->ReleaseLongArrayElements(destination, host, JNI_ABORT);
        cudaFree(output);
        jclass cls = env->FindClass("java/lang/IllegalStateException");
        if (cls) env->ThrowNew(cls, error.what());
    }
}

extern "C" JNIEXPORT void JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_fillFocusOn(JNIEnv* env, jclass,
    jlong handle, jint x, jintArray jSelected, jlongArray destination) {
    auto* index = reinterpret_cast<FocusIndex*>(handle);
    int* dMap = nullptr;
    unsigned long long* output = nullptr;
    auto release = [&] { cudaFree(dMap); cudaFree(output); };
    try {
        if (!index || x < 0 || x >= index->n || jSelected == nullptr)
            throw std::invalid_argument("invalid S6 backbone focus request");
        int m = env->GetArrayLength(jSelected);
        if (m < 2 || m >= index->n ||
                env->GetArrayLength(destination) != 2LL * m * m)
            throw std::invalid_argument("invalid S6 backbone focus dimensions");
        std::vector<jint> selected(m);
        env->GetIntArrayRegion(jSelected, 0, m, selected.data());
        if (env->ExceptionCheck()) return;
        std::vector<int> map(index->n, -1);
        for (int i = 0; i < m; i++) {
            int taxon = selected[i];
            if (taxon < 0 || taxon >= index->n || taxon == x || map[taxon] >= 0)
                throw std::invalid_argument("invalid S6 backbone taxon");
            map[taxon] = i;
        }
        if (index->staged) {
            std::vector<jlong> host(2ULL * m * m);
            SelectedGuard guard(index->selectedBusy);
            prepareSelected(*index, selected, map);
            runCompactFocus(*index, x, index->selectedChildren,
                index->selectedTaxa, map, true, m,
                reinterpret_cast<unsigned long long*>(host.data()));
            env->SetLongArrayRegion(destination, 0,
                static_cast<jsize>(host.size()), host.data());
            return;
        }
        size_t bytes = 2ULL * m * m * sizeof(uint64_t);
        check(cudaMalloc(&dMap, index->n * sizeof(int)),
            "allocating backbone taxon map");
        check(cudaMalloc(&output, bytes), "allocating backbone focus");
        check(cudaMemcpy(dMap, map.data(), index->n * sizeof(int),
            cudaMemcpyHostToDevice), "copying backbone taxon map");
        check(cudaMemset(output, 0, bytes), "clearing backbone focus");
        int first = index->offsets[x], last = index->offsets[x + 1];
        if (last > first) {
            int blocks = (last - first + 255) / 256;
            focusedSubsetByMembership<<<blocks, 256>>>(m,
                index->parts, index->children, index->taxa,
                index->partOfChild, index->taxonChildren, dMap,
                first, last, output);
            check(cudaGetLastError(), "launching backbone focus kernel");
        }
        std::vector<jlong> host(2ULL * m * m);
        check(cudaMemcpy(host.data(), output, bytes, cudaMemcpyDeviceToHost),
            "copying backbone focus");
        env->SetLongArrayRegion(destination, 0, static_cast<jsize>(host.size()),
            host.data());
        release();
    } catch (const std::exception& error) {
        release();
        jclass cls = env->FindClass("java/lang/IllegalStateException");
        if (cls) env->ThrowNew(cls, error.what());
    }
}

extern "C" JNIEXPORT void JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_destroyFocusIndex(JNIEnv*, jclass,
    jlong handle) {
    delete reinterpret_cast<FocusIndex*>(handle);
}
