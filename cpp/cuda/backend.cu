// Standalone CUDA backend. All ownership stays in C++; no JNI arrays or objects
// cross the device boundary. Original species sets are packed, duplicate-free
// bitsets, so multicopy leaves never multiply an intersection or a triplet.
#include "stelar/core.hpp"
#include <climits>
#include <cmath>
#include <cstring>
#include <cuda_runtime.h>
#include <numeric>
#include <sstream>

namespace stelar {
namespace {
using U64 = unsigned long long;
constexpr unsigned block_size = 256;
struct U128 {
    U64 lo, hi;
};
static_assert(sizeof(U128) == 16, "two-limb score layout");

void check(cudaError_t e, const char *what) {
    if (e != cudaSuccess)
        throw std::runtime_error(std::string("CUDA ") + what + ": " + cudaGetErrorString(e));
}
size_t multiply_size(size_t a, size_t b) {
    if (b && a > std::numeric_limits<size_t>::max() / b)
        throw std::length_error("CUDA array size overflow");
    return a * b;
}
template <class T> class DeviceArray {
    T *ptr_ = nullptr;
    size_t capacity_ = 0;

  public:
    DeviceArray() = default;
    explicit DeviceArray(size_t n) {
        reserve(n);
    }
    ~DeviceArray() {
        if (ptr_)
            cudaFree(ptr_);
    }
    DeviceArray(const DeviceArray &) = delete;
    DeviceArray &operator=(const DeviceArray &) = delete;
    T *data() const {
        return ptr_;
    }
    size_t bytes() const {
        return capacity_ * sizeof(T);
    }
    void reserve(size_t n) {
        if (n <= capacity_)
            return;
        // These buffers are scratch or are overwritten by upload; preserving
        // their discarded contents would briefly double the live device RAM.
        if (ptr_)
            check(cudaFree(ptr_), "release resized device array");
        ptr_ = nullptr;
        capacity_ = 0;
        check(cudaMalloc(reinterpret_cast<void **>(&ptr_), multiply_size(n, sizeof(T))),
              "allocate device array");
        capacity_ = n;
    }
    void upload(const T *p, size_t n) {
        reserve(n);
        if (n)
            check(cudaMemcpy(ptr_, p, n * sizeof(T), cudaMemcpyHostToDevice),
                  "upload device array");
    }
    void upload(const std::vector<T> &v) {
        upload(v.data(), v.size());
    }
    void download(T *p, size_t n) const {
        if (n)
            check(cudaMemcpy(p, ptr_, n * sizeof(T), cudaMemcpyDeviceToHost),
                  "download device array");
    }
    void clear(size_t n) {
        reserve(n);
        if (n)
            check(cudaMemset(ptr_, 0, n * sizeof(T)), "clear device array");
    }
};
size_t budget(const Config &cfg, int cap_mb = 0) {
    size_t free = 0, total = 0;
    check(cudaMemGetInfo(&free, &total), "query free device memory");
    const double fraction = std::min(.98, std::max(.01, cfg.gpu_occupancy));
    size_t bytes = static_cast<size_t>(free * fraction);
    if (cap_mb > 0)
        bytes = std::min(bytes, multiply_size(size_t(cap_mb), size_t(1) << 20));
    if (bytes < (size_t(1) << 20))
        throw std::runtime_error("less than 1 MiB available within CUDA memory budget");
    return bytes;
}
using GSet = SetPool::Record;
using GPart = Observation;
struct PackedObservations {
    std::vector<GSet> sets;
    std::vector<U64> bits;
    std::vector<Id> children;
    std::vector<GPart> parts;
};
PackedObservations pack_observations(const Observations &obs, size_t begin, size_t end,
                                     bool whole_pool) {
    PackedObservations p;
    if (whole_pool) {
        p.sets.reserve(obs.sets.records.size());
        for (const auto &r : obs.sets.records)
            p.sets.push_back({r.offset, r.words});
        p.bits.assign(obs.sets.arena.begin(), obs.sets.arena.end());
        p.children = obs.child_sets;
        p.parts.reserve(end - begin);
        for (size_t i = begin; i < end; ++i) {
            const auto &v = obs.partitions[i];
            p.parts.push_back({v.offset, v.children, v.frequency});
        }
        return p;
    }
    std::unordered_map<Id, Id> mapping;
    for (size_t i = begin; i < end; ++i) {
        const auto &v = obs.partitions[i];
        p.parts.push_back({p.children.size(), v.children, v.frequency});
        for (uint32_t j = 0; j < v.children; ++j) {
            Id old = obs.child_sets[v.offset + j];
            auto it = mapping.find(old);
            Id id;
            if (it == mapping.end()) {
                if (p.sets.size() >= none)
                    throw std::length_error("too many CUDA sets");
                id = static_cast<Id>(p.sets.size());
                mapping.emplace(old, id);
                const auto &r = obs.sets.records[old];
                p.sets.push_back({p.bits.size(), r.words, r.count});
                const auto *data = obs.sets.data(old);
                p.bits.insert(p.bits.end(), data, data + r.words);
            } else
                id = it->second;
            p.children.push_back(id);
        }
    }
    return p;
}
size_t resident_bytes(const Observations &obs) {
    return obs.sets.records.size() * sizeof(GSet) + obs.sets.arena.size() * sizeof(U64) +
           obs.child_sets.size() * sizeof(Id) + obs.partitions.size() * sizeof(GPart);
}
struct DeviceObservations {
    DeviceArray<GSet> sets;
    DeviceArray<U64> bits;
    DeviceArray<Id> children;
    DeviceArray<GPart> parts;
    void upload(const PackedObservations &p) {
        sets.upload(p.sets);
        bits.upload(p.bits);
        children.upload(p.children);
        parts.upload(p.parts);
    }
    void upload_whole(const Observations &p) {
        // The C++ host arena is already the device's POD representation.
        // cudaMemcpy reads bytes; no second host forest/bitset arena is built.
        sets.upload(p.sets.records);
        bits.upload(reinterpret_cast<const U64 *>(p.sets.arena.data()), p.sets.arena.size());
        children.upload(p.child_sets);
        parts.upload(p.partitions);
    }
};
// Conservative ranges bound every intermediate doubled weight before launch.
// Division first avoids overflow in the bound itself.
Score maximum_two_score(const Observations &obs, size_t n) {
    Score freq = 0;
    for (const auto &p : obs.partitions)
        freq += p.frequency;
    if (n < 3)
        return 0;
    uint64_t a = n, b = n - 1, c = n - 2;
    if (a % 3 == 0)
        a /= 3;
    else if (b % 3 == 0)
        b /= 3;
    else
        c /= 3;
    Score answer = freq, max = ~Score(0);
    for (uint64_t x : {a, b, c}) {
        if (x && answer > max / x)
            throw std::overflow_error("triplet score may exceed exact 128-bit range");
        answer *= x;
    }
    return answer;
}
__device__ __forceinline__ U128 plus128(U128 a, U128 b) {
    U64 lo = a.lo + b.lo;
    return {lo, a.hi + b.hi + (lo < a.lo)};
}
__device__ __forceinline__ U128 product128(U64 a, U64 b) {
    return {a * b, __umul64hi(a, b)};
}
__device__ __forceinline__ U128 scale128(U128 a, U64 f) {
    U128 b = product128(a.lo, f);
    b.hi += a.hi * f;
    return b;
}
__device__ __forceinline__ uint32_t intersection(const U64 *a, uint32_t aw, const U64 *bits,
                                                 GSet s) {
    uint32_t result = 0, len = min(aw, s.words);
    for (uint32_t w = 0; w < len; ++w)
        result += __popcll(a[w] & bits[s.offset + w]);
    return result;
}
__device__ __forceinline__ bool has(const U64 *bits, GSet s, Id x) {
    uint32_t w = x >> 6;
    return w < s.words && ((bits[s.offset + w] >> (x & 63)) & 1ULL);
}
template <bool Wide> struct Arithmetic;
template <> struct Arithmetic<true> {
    using Acc = U128;
    __device__ static Acc zero() {
        return {0, 0};
    }
    __device__ static Acc add(Acc a, Acc b) {
        return plus128(a, b);
    }
    __device__ static Acc triple(U64 a, U64 b, U64 c) {
        return scale128(product128(a, b), c);
    }
    __device__ static Acc scale(Acc a, U64 b) {
        return scale128(a, b);
    }
};
template <> struct Arithmetic<false> {
    using Acc = U64;
    __device__ static Acc zero() {
        return 0;
    }
    __device__ static Acc add(Acc a, Acc b) {
        return a + b;
    }
    __device__ static Acc triple(U64 a, U64 b, U64 c) {
        return a * b * c;
    }
    __device__ static Acc scale(Acc a, U64 b) {
        return a * b;
    }
};
template <bool Wide>
__global__ void weight_kernel(const U64 *candidate, uint32_t words, size_t split_count,
                              const GSet *sets, const U64 *bits, const Id *children,
                              const GPart *parts, size_t part_count,
                              typename Arithmetic<Wide>::Acc *output) {
    using Op = Arithmetic<Wide>;
    using Acc = typename Op::Acc;
    __shared__ Acc sums[block_size];
    size_t split = blockIdx.x;
    if (split >= split_count)
        return;
    const U64 *a = candidate + (2 * split) * words;
    const U64 *b = a + words;
    Acc sum = Op::zero();
    for (size_t i = threadIdx.x; i < part_count; i += blockDim.x) {
        GPart p = parts[i];
        Acc value = Op::zero();
        if (p.children == 2) {
            GSet c0 = sets[children[p.offset]], c1 = sets[children[p.offset + 1]];
            U64 a0 = intersection(a, words, bits, c0), a1 = intersection(a, words, bits, c1);
            U64 b0 = intersection(b, words, bits, c0), b1 = intersection(b, words, bits, c1);
            if (a0 && b1)
                value = Op::add(value, Op::triple(a0, b1, a0 + b1 - 2));
            if (a1 && b0)
                value = Op::add(value, Op::triple(a1, b0, a1 + b0 - 2));
        } else {
            U64 sa = 0, sb = 0;
            for (uint32_t j = 0; j < p.children; ++j) {
                GSet s = sets[children[p.offset + j]];
                sa += intersection(a, words, bits, s);
                sb += intersection(b, words, bits, s);
            }
            for (uint32_t j = 0; j < p.children; ++j) {
                GSet s = sets[children[p.offset + j]];
                U64 ai = intersection(a, words, bits, s), bi = intersection(b, words, bits, s);
                if (ai >= 2 && sb > bi)
                    value = Op::add(value, Op::triple(ai, ai - 1, sb - bi));
                if (bi >= 2 && sa > ai)
                    value = Op::add(value, Op::triple(bi, bi - 1, sa - ai));
            }
        }
        sum = Op::add(sum, Op::scale(value, p.frequency));
    }
    sums[threadIdx.x] = sum;
    __syncthreads();
    for (unsigned stride = blockDim.x / 2; stride; stride >>= 1) {
        if (threadIdx.x < stride)
            sums[threadIdx.x] = Op::add(sums[threadIdx.x], sums[threadIdx.x + stride]);
        __syncthreads();
    }
    if (!threadIdx.x)
        output[split] = sums[0];
}
size_t choose_split_batch(size_t q, size_t per_split, size_t available, size_t resident,
                          const Config &cfg) {
    size_t requested = q;
    if (!cfg.no_gpu_batch) {
        if (cfg.gpu_batches > 0)
            requested = (q + cfg.gpu_batches - 1) / cfg.gpu_batches;
        else if (cfg.gpu_batch_size > 0)
            requested = cfg.gpu_batch_size;
        else if (cfg.gpu_control > 0)
            requested = std::max(size_t(1), size_t(resident * cfg.gpu_control) / per_split);
        else
            requested = std::max(size_t(1), available / per_split);
    }
    size_t fit = available / per_split;
    if (!fit)
        throw std::runtime_error("CUDA memory budget cannot hold a single candidate split");
    if ((cfg.no_gpu_batch || cfg.gpu_batches > 0 || cfg.gpu_batch_size > 0) && requested > fit)
        throw std::runtime_error(
            "explicit GPU split batch exceeds available CUDA memory; reduce batch size");
    return std::min({q, requested, fit, size_t(INT_MAX)});
}
void pack_candidates(const CandidateSpace &space, const std::vector<Id> &ids, size_t begin,
                     size_t count, size_t words, std::vector<U64> &packed) {
    packed.assign(multiply_size(count, 2 * words), 0);
    for (size_t i = 0; i < count; ++i) {
        const auto &s = space.splits.at(ids[begin + i]);
        Id sides[2] = {s.left, s.right};
        for (size_t j = 0; j < 2; ++j) {
            const auto &r = space.sets.records.at(sides[j]);
            std::copy_n(space.sets.data(sides[j]), r.words, packed.data() + (2 * i + j) * words);
        }
    }
}

// Leaf-only pair queries use an Euler RMQ. The low 32 bits carry subtree size;
// the high 32 bits carry depth, making unsigned minimum return the LCA size.
struct SimTree {
    uint64_t leaf_offset, rmq_offset;
    uint32_t leaves, euler;
};
struct SimBatch {
    std::vector<SimTree> trees;
    std::vector<Id> taxa;
    std::vector<uint32_t> first;
    std::vector<U64> rmq;
    size_t bytes() const {
        return trees.size() * sizeof(SimTree) + taxa.size() * sizeof(Id) +
               first.size() * sizeof(uint32_t) + rmq.size() * sizeof(U64);
    }
};
void append_sim_tree(const Tree &t, SimBatch &batch) {
    auto post = t.postorder();
    if (t.root == none || t.leaves.size() < 3)
        return;
    std::vector<uint32_t> sizes(t.nodes.size());
    for (Id v : post) {
        if (t.nodes[v].taxon != none)
            sizes[v] = 1;
        else
            for (Id c = t.nodes[v].child; c != none; c = t.nodes[c].next)
                sizes[v] += sizes[c];
    }
    const uint32_t k = sizes[t.root];
    if (k < 3)
        return;
    std::vector<U64> euler;
    std::vector<std::pair<Id, uint32_t>> leaves;
    struct Frame {
        Id node, child;
        uint32_t depth;
    };
    std::vector<Frame> stack;
    auto enter = [&](Id v, uint32_t depth) {
        if (t.nodes[v].taxon != none)
            leaves.emplace_back(t.nodes[v].taxon, static_cast<uint32_t>(euler.size()));
        euler.push_back((U64(depth) << 32) | sizes[v]);
        stack.push_back({v, t.nodes[v].child, depth});
    };
    enter(t.root, 0);
    while (!stack.empty()) {
        Frame &f = stack.back();
        if (f.child != none) {
            Id c = f.child;
            f.child = t.nodes[c].next;
            uint32_t d = f.depth + 1;
            enter(c, d);
        } else {
            stack.pop_back();
            if (!stack.empty()) {
                const Frame &p = stack.back();
                euler.push_back((U64(p.depth) << 32) | sizes[p.node]);
            }
        }
    }
    if (euler.size() > UINT32_MAX)
        throw std::length_error("single CUDA similarity tree too large");
    std::sort(leaves.begin(), leaves.end());
    for (size_t i = 1; i < leaves.size(); ++i)
        if (leaves[i - 1].first == leaves[i].first)
            throw std::invalid_argument("triplet similarity expects single-copy DISCO trees");
    SimTree desc{batch.taxa.size(), batch.rmq.size(), k, static_cast<uint32_t>(euler.size())};
    for (const auto &l : leaves) {
        batch.taxa.push_back(l.first);
        batch.first.push_back(l.second);
    }
    batch.rmq.insert(batch.rmq.end(), euler.begin(), euler.end());
    size_t length = euler.size();
    for (size_t span = 2; span <= length; span *= 2) {
        size_t previous = batch.rmq.size() - length;
        for (size_t i = 0; i < length; ++i)
            batch.rmq.push_back(std::min(batch.rmq[previous + i],
                                         batch.rmq[previous + std::min(length - 1, i + span / 2)]));
        if (span > length / 2)
            break;
    }
    batch.trees.push_back(desc);
}
__device__ int leaf_position(Id x, SimTree t, const Id *taxa, const uint32_t *first) {
    uint32_t lo = 0, hi = t.leaves;
    while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (taxa[t.leaf_offset + mid] < x)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo < t.leaves && taxa[t.leaf_offset + lo] == x
               ? static_cast<int>(first[t.leaf_offset + lo])
               : -1;
}
__device__ void decode_pair(uint64_t p, uint32_t n, uint32_t &a, uint32_t &b) {
    // Match Matrix's row-major upper triangle exactly. Integer binary search
    // also avoids rounding near row boundaries for very large species counts.
    uint64_t lo = 0, hi = n - 1;
    while (lo + 1 < hi) {
        uint64_t mid = lo + (hi - lo) / 2, start = mid * n - mid * (mid + 1) / 2;
        if (start <= p)
            lo = mid;
        else
            hi = mid;
    }
    a = static_cast<uint32_t>(lo);
    b = static_cast<uint32_t>(lo + 1 + p - (lo * n - lo * (lo + 1) / 2));
}
__global__ void similarity_kernel(uint64_t start, size_t count, uint32_t n, const SimTree *trees,
                                  size_t tree_count, const Id *taxa, const uint32_t *first,
                                  const U64 *rmq, U64 *numerator, U64 *denominator) {
    size_t at = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (at >= count)
        return;
    uint32_t a, b;
    decode_pair(start + at, n, a, b);
    U64 num = 0, den = 0;
    for (size_t i = 0; i < tree_count; ++i) {
        SimTree t = trees[i];
        int p = leaf_position(a, t, taxa, first), q = leaf_position(b, t, taxa, first);
        if (p < 0 || q < 0)
            continue;
        if (p > q) {
            int tmp = p;
            p = q;
            q = tmp;
        }
        uint32_t length = q - p + 1, level = 31 - __clz(length), span = 1U << level;
        uint64_t off = t.rmq_offset + uint64_t(level) * t.euler;
        U64 lca = min(rmq[off + p], rmq[off + q - span + 1]);
        num += t.leaves - static_cast<uint32_t>(lca);
        den += t.leaves - 2;
    }
    numerator[at] = num;
    denominator[at] = den;
}

// A selected-set cache contains O(S ceil(m/64)) bits, not an n^3 tensor.
// The outgroup row convention mirrors the Java focus accessor exactly.
__global__ void select_sets_kernel(const GSet *sets, const U64 *bits, size_t count,
                                   const Id *selected, uint32_t m, uint32_t words, U64 *output) {
    size_t at = size_t(blockIdx.x) * blockDim.x + threadIdx.x, total = count * words;
    if (at >= total)
        return;
    size_t s = at / words;
    uint32_t w = at % words;
    U64 value = 0;
    for (uint32_t j = w * 64; j < min(m, (w + 1) * 64); ++j)
        if (has(bits, sets[s], selected[j]))
            value |= 1ULL << (j & 63);
    output[at] = value;
}
__device__ void add_focus_partition(const GPart p, const Id *children, const U64 *selected,
                                    const Id *selected_ids, uint32_t words, uint32_t m,
                                    uint32_t side, U64 *pair, U64 *sister) {
    const U64 *xb = selected + size_t(children[p.offset + side]) * words;
    for (uint32_t out = 0; out < p.children; ++out) {
        if (out == side)
            continue;
        const U64 *ob = selected + size_t(children[p.offset + out]) * words;
        for (uint32_t wa = 0; wa < words; ++wa) {
            U64 aa = ob[wa];
            while (aa) {
                unsigned bit = __ffsll(aa) - 1;
                aa &= aa - 1;
                uint32_t a = wa * 64 + bit;
                // Enumerate once and preserve Java's GLOBAL-ID upper orientation,
                // even when the selected backbone is in an unsorted DFS order.
                for (uint32_t wb = wa; wb < words; ++wb) {
                    U64 bb = ob[wb];
                    if (wb == wa)
                        bb &= (bit == 63 ? 0ULL : (~0ULL << (bit + 1)));
                    while (bb) {
                        unsigned bbit = __ffsll(bb) - 1;
                        bb &= bb - 1;
                        uint32_t b = wb * 64 + bbit;
                        size_t cell = selected_ids[a] < selected_ids[b] ? size_t(a) * m + b
                                                                        : size_t(b) * m + a;
                        atomicAdd(pair + cell, U64(p.frequency));
                    }
                }
                for (uint32_t wb = 0; wb < words; ++wb) {
                    U64 bb = xb[wb];
                    while (bb) {
                        unsigned bbit = __ffsll(bb) - 1;
                        bb &= bb - 1;
                        uint32_t b = wb * 64 + bbit;
                        if (a != b)
                            atomicAdd(sister + size_t(a) * m + b, U64(p.frequency));
                    }
                }
            }
        }
    }
}
__global__ void focus_kernel(const GSet *sets, const U64 *bits, const Id *children,
                             const GPart *parts, size_t count, const U64 *selected,
                             const Id *selected_ids, uint32_t words, uint32_t m, Id x,
                             U64 *output) {
    size_t at = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (at >= count)
        return;
    GPart p = parts[at];
    for (uint32_t side = 0; side < p.children; ++side)
        if (has(bits, sets[children[p.offset + side]], x))
            add_focus_partition(p, children, selected, selected_ids, words, m, side, output,
                                output + size_t(m) * m);
}
} // namespace

bool cuda_available(std::string *reason) {
    int count = 0;
    cudaError_t e = cudaGetDeviceCount(&count);
    if (e != cudaSuccess || count == 0) {
        if (reason)
            *reason = e == cudaSuccess ? "no CUDA device" : cudaGetErrorString(e);
        cudaGetLastError();
        return false;
    }
    if (reason) {
        cudaDeviceProp properties{};
        int device = 0;
        if (cudaGetDevice(&device) == cudaSuccess &&
            cudaGetDeviceProperties(&properties, device) == cudaSuccess)
            *reason = std::string(properties.name) + " (compute " +
                      std::to_string(properties.major) + "." + std::to_string(properties.minor) +
                      ")";
        else
            *reason = "CUDA device available";
    }
    return true;
}

std::vector<Score> cuda_weights(const CandidateSpace &space, const Observations &obs,
                                const std::vector<Id> &ids, const Config &cfg) {
    std::vector<Score> result(ids.size(), 0);
    if (ids.empty() || obs.partitions.empty())
        return result;
    Score bound = maximum_two_score(obs, space.taxa);
    bool wide = bound > Score(LLONG_MAX) || cfg.score_type == "int128" ||
                cfg.score_type == "i128" || cfg.score_type == "exact";
    const size_t words = (space.taxa + 63) / 64, per_split = 2 * words * sizeof(U64) + sizeof(U128);
    const size_t available = budget(cfg), resident = resident_bytes(obs);
    bool full = resident < available * 2 / 3;
    // Oversized observations are staged by whole partitions. The only scores
    // retained across stages are one exact accumulator per candidate split.
    std::vector<std::pair<size_t, size_t>> stages;
    if (full)
        stages.emplace_back(0, obs.partitions.size());
    else {
        size_t start = 0, bytes = 0, limit = available / 2;
        for (size_t i = 0; i < obs.partitions.size(); ++i) {
            const auto &p = obs.partitions[i];
            size_t need = sizeof(GPart) + p.children * (sizeof(Id) + sizeof(GSet));
            for (uint32_t j = 0; j < p.children; ++j)
                need += obs.sets.records[obs.child_sets[p.offset + j]].words * sizeof(U64);
            if (need > limit)
                throw std::runtime_error("one observation exceeds CUDA staging budget");
            if (bytes && need > limit - bytes) {
                stages.emplace_back(start, i);
                start = i;
                bytes = 0;
            }
            bytes += need;
        }
        stages.emplace_back(start, obs.partitions.size());
    }
    log(cfg, "CUDA exact triplet weights: " + std::to_string(ids.size()) + " splits, " +
                 std::to_string(obs.partitions.size()) + " partitions, " +
                 (wide ? "128-bit" : "64-bit") + ", " + std::to_string(stages.size()) +
                 " observation stage(s)");
    DeviceObservations device;
    DeviceArray<U64> candidates;
    DeviceArray<U128> scores128;
    DeviceArray<U64> scores64;
    std::vector<U64> packed, host64;
    std::vector<U128> host128;
    auto last = std::chrono::steady_clock::now();
    size_t completed = 0, total = stages.size() * ids.size();
    for (const auto &stage : stages) {
        PackedObservations p;
        if (full)
            device.upload_whole(obs);
        else {
            p = pack_observations(obs, stage.first, stage.second, false);
            device.upload(p);
        }
        const size_t partition_count = stage.second - stage.first;
        size_t current = budget(cfg) + candidates.bytes() + scores128.bytes() + scores64.bytes();
        size_t batch = choose_split_batch(ids.size(), per_split, current, resident, cfg);
        candidates.reserve(2 * words * batch);
        if (wide)
            scores128.reserve(batch);
        else
            scores64.reserve(batch);
        for (size_t begin = 0; begin < ids.size(); begin += batch) {
            size_t count = std::min(batch, ids.size() - begin);
            pack_candidates(space, ids, begin, count, words, packed);
            candidates.upload(packed);
            if (wide) {
                weight_kernel<true><<<static_cast<unsigned>(count), block_size>>>(
                    candidates.data(), words, count, device.sets.data(), device.bits.data(),
                    device.children.data(), device.parts.data(), partition_count, scores128.data());
                check(cudaGetLastError(), "launch exact 128-bit weights");
                host128.resize(count);
                scores128.download(host128.data(), count);
                for (size_t i = 0; i < count; ++i)
                    result[begin + i] += (Score(host128[i].hi) << 64) | host128[i].lo;
            } else {
                weight_kernel<false><<<static_cast<unsigned>(count), block_size>>>(
                    candidates.data(), words, count, device.sets.data(), device.bits.data(),
                    device.children.data(), device.parts.data(), partition_count, scores64.data());
                check(cudaGetLastError(), "launch exact 64-bit weights");
                host64.resize(count);
                scores64.download(host64.data(), count);
                for (size_t i = 0; i < count; ++i)
                    result[begin + i] += host64[i];
            }
            completed += count;
            auto now = std::chrono::steady_clock::now();
            if (cfg.gpu_progress > 0 &&
                std::chrono::duration<double>(now - last).count() >= cfg.gpu_progress) {
                log(cfg, "CUDA weight progress " + std::to_string(completed) + "/" +
                             std::to_string(total));
                last = now;
            }
        }
    }
    for (auto &score : result) {
        if (score & 1)
            throw std::runtime_error("nonintegral doubled CUDA triplet weight");
        score >>= 1;
    }
    return result;
}

bool cuda_triplet_matrix(const std::vector<Tree> &trees, size_t n, const Config &cfg,
                         Matrix &output) {
    if (!cuda_available())
        return false;
    if (n > UINT32_MAX)
        throw std::length_error("too many CUDA taxa");
    output = Matrix(n);
    if (n < 2 || trees.empty())
        return true;
    size_t limit = budget(cfg, cfg.gpu_sim_cap_mb), pairs = output.values.size();
    size_t tile = std::min(pairs, std::max(size_t(1), (limit / 3) / (2 * sizeof(U64))));
    // gpu_tile is the edge length of the pair output tile; zero selects budget.
    if (cfg.gpu_tile > 0)
        tile = std::min(tile, multiply_size(size_t(cfg.gpu_tile), size_t(cfg.gpu_tile)));
    DeviceArray<U64> num(tile), den(tile), rmq;
    DeviceArray<SimTree> desc;
    DeviceArray<Id> taxa;
    DeviceArray<uint32_t> first;
    std::vector<U64> hnum(tile), hden(tile), denominators(pairs, 0);
    auto consume = [&](const SimBatch &batch) {
        if (batch.trees.empty())
            return;
        desc.upload(batch.trees);
        taxa.upload(batch.taxa);
        first.upload(batch.first);
        rmq.upload(batch.rmq);
        for (size_t begin = 0; begin < pairs; begin += tile) {
            size_t count = std::min(tile, pairs - begin);
            similarity_kernel<<<static_cast<unsigned>((count + block_size - 1) / block_size),
                                block_size>>>(begin, count, n, desc.data(), batch.trees.size(),
                                              taxa.data(), first.data(), rmq.data(), num.data(),
                                              den.data());
            check(cudaGetLastError(), "launch triplet similarity");
            num.download(hnum.data(), count);
            den.download(hden.data(), count);
            for (size_t i = 0; i < count; ++i) {
                output.values[begin + i] += double(hnum[i]);
                denominators[begin + i] += hden[i];
            }
        }
    };
    SimBatch batch;
    for (const auto &tree : trees) {
        size_t k = tree.leaves.size();
        if (k < 3)
            continue;
        size_t length = tree.nodes.size() * 2, levels = 1;
        for (size_t v = length; v > 1; v >>= 1)
            ++levels;
        size_t estimate = sizeof(SimTree) + k * (sizeof(Id) + sizeof(uint32_t)) +
                          multiply_size(length, levels * sizeof(U64));
        if (estimate > limit * 2 / 3)
            throw std::runtime_error("one similarity tree exceeds CUDA staging cap; use CPU or "
                                     "raise --gpu-sim-vram-cap-mb");
        if (!batch.trees.empty() && batch.bytes() + estimate > limit * 2 / 3) {
            consume(batch);
            batch = SimBatch();
        }
        append_sim_tree(tree, batch);
    }
    consume(batch);
    for (size_t i = 0; i < pairs; ++i)
        output.values[i] = denominators[i] ? output.values[i] / double(denominators[i]) : 0;
    return true;
}

namespace {
class CudaFocusImpl final : public CudaFocus {
    const Observations &observations_;
    size_t n_;
    Config config_;
    bool resident_;
    DeviceObservations device_;
    DeviceArray<Id> selected_ids_;
    DeviceArray<U64> selected_bits_, output_;
    std::vector<Id> selected_cache_;
    std::vector<U64> host_;
    std::mutex mutex_;
    void prepare_selected(const std::vector<Id> &selected, size_t set_count, uint32_t words) {
        selected_ids_.upload(selected);
        selected_bits_.reserve(multiply_size(set_count, words));
        size_t total = multiply_size(set_count, words);
        if (total) {
            select_sets_kernel<<<static_cast<unsigned>((total + block_size - 1) / block_size),
                                 block_size>>>(device_.sets.data(), device_.bits.data(), set_count,
                                               selected_ids_.data(), selected.size(), words,
                                               selected_bits_.data());
            check(cudaGetLastError(), "prepare selected species bitsets");
        }
    }
    void launch(size_t partitions, Id x, uint32_t words, uint32_t m) {
        if (!partitions)
            return;
        focus_kernel<<<static_cast<unsigned>((partitions + block_size - 1) / block_size),
                       block_size>>>(
            device_.sets.data(), device_.bits.data(), device_.children.data(), device_.parts.data(),
            partitions, selected_bits_.data(), selected_ids_.data(), words, m, x, output_.data());
        check(cudaGetLastError(), "launch focused triplet support");
    }

  public:
    CudaFocusImpl(const Observations &o, size_t n, const Config &c)
        : observations_(o), n_(n), config_(c), resident_(false) {
        // Keep space for a selected backbone and focused output when deciding
        // whether the original packed observations can remain on the device.
        resident_ = resident_bytes(o) < budget(c, c.gpu_treewalk_cap_mb) / 2;
        if (resident_)
            device_.upload_whole(o);
    }
    Focus focus(Id x, const std::vector<Id> &selected) override {
        std::lock_guard<std::mutex> guard(mutex_);
        if (x >= n_)
            throw std::out_of_range("focused taxon outside registry");
        Focus answer;
        answer.n = n_;
        answer.selected = selected;
        size_t m = selected.size();
        size_t cells = multiply_size(m, m);
        answer.pair.resize(cells);
        answer.sister.resize(cells);
        if (!m)
            return answer;
        std::vector<uint8_t> seen(n_);
        for (Id id : selected) {
            if (id >= n_ || id == x || seen[id])
                throw std::invalid_argument("invalid or repeated selected focus species");
            seen[id] = 1;
        }
        size_t bytes = multiply_size(cells, 2 * sizeof(U64)),
               allowed = budget(config_, config_.gpu_treewalk_cap_mb);
        if (bytes > allowed / 2)
            throw std::runtime_error(
                "focused support exceeds CUDA cap; reduce backbone size or raise cap");
        uint32_t words = (m + 63) / 64;
        output_.clear(cells * 2);
        if (resident_) {
            if (selected_cache_ != selected) {
                size_t cache =
                    multiply_size(observations_.sets.records.size(), words * sizeof(U64));
                if (cache > allowed / 2)
                    throw std::runtime_error("selected focus cache exceeds CUDA cap");
                prepare_selected(selected, observations_.sets.records.size(), words);
                selected_cache_ = selected;
            }
            launch(observations_.partitions.size(), x, words, m);
        } else {
            // Compact only relevant original partitions, in bounded chunks.
            // The species sets preserve every multiplicity through frequency.
            const size_t stage_limit = allowed / 2;
            size_t begin = 0, used = 0;
            auto consume = [&](size_t lo, size_t hi) {
                if (lo == hi)
                    return;
                Observations filtered;
                for (size_t i = lo; i < hi; ++i) {
                    const auto &part = observations_.partitions[i];
                    bool relevant = false;
                    for (uint32_t j = 0; j < part.children && !relevant; ++j) {
                        Id s = observations_.child_sets[part.offset + j];
                        const auto &r = observations_.sets.records[s];
                        relevant = (x >> 6) < r.words &&
                                   ((observations_.sets.data(s)[x >> 6] >> (x & 63)) & 1);
                    }
                    if (!relevant)
                        continue;
                    Observation p{filtered.child_sets.size(), part.children, part.frequency};
                    for (uint32_t j = 0; j < part.children; ++j) {
                        Id s = observations_.child_sets[part.offset + j];
                        filtered.child_sets.push_back(
                            filtered.sets.intern(observations_.sets.bits(s)));
                    }
                    filtered.partitions.push_back(p);
                }
                if (filtered.partitions.empty())
                    return;
                device_.upload_whole(filtered);
                prepare_selected(selected, filtered.sets.records.size(), words);
                launch(filtered.partitions.size(), x, words, m);
            };
            for (size_t i = 0; i < observations_.partitions.size(); ++i) {
                const auto &p = observations_.partitions[i];
                size_t need =
                    sizeof(GPart) + p.children * (sizeof(Id) + sizeof(GSet) + words * sizeof(U64));
                for (uint32_t j = 0; j < p.children; ++j)
                    need +=
                        observations_.sets.records[observations_.child_sets[p.offset + j]].words *
                        sizeof(U64);
                if (need > stage_limit)
                    throw std::runtime_error("one focused partition exceeds CUDA staging budget");
                if (used && need > stage_limit - used) {
                    consume(begin, i);
                    begin = i;
                    used = 0;
                }
                used += need;
            }
            consume(begin, observations_.partitions.size());
        }
        host_.resize(cells * 2);
        output_.download(host_.data(), host_.size());
        std::copy_n(host_.data(), cells, answer.pair.data());
        std::copy_n(host_.data() + cells, cells, answer.sister.data());
        return answer;
    }
};
} // namespace
std::unique_ptr<CudaFocus> make_cuda_focus(const Observations &obs, size_t n, const Config &cfg) {
    if (!cuda_available())
        return nullptr;
    maximum_two_score(obs, n);
    return std::make_unique<CudaFocusImpl>(obs, n, cfg);
}

namespace {
struct CrossHash {
    U64 sum;
    Id id;
    uint32_t count;
};
struct CrossResult {
    Id parent, left, right;
};
__global__ void cross_kernel(const GSet *sets, const U64 *words, const CrossHash *hashes,
                             size_t hash_count, const Id *parents, size_t parent_begin,
                             size_t parent_count, const Id *sides, const uint32_t *counts,
                             const U64 *sums, size_t side_begin, size_t side_end,
                             CrossResult *result, U64 *emitted) {
    size_t pi = blockIdx.x;
    if (pi >= parent_count)
        return;
    Id parent = parents[parent_begin + pi];
    GSet a = sets[parent];
    for (size_t k = side_begin + threadIdx.x; k < side_end; k += blockDim.x) {
        Id left = sides[k];
        if (counts[left] > counts[parent] / 2)
            break;
        GSet b = sets[left];
        bool contained = true;
        for (uint32_t w = 0; w < b.words; ++w)
            if (words[b.offset + w] & ~(w < a.words ? words[a.offset + w] : 0ULL)) {
                contained = false;
                break;
            }
        if (!contained)
            continue;
        U64 target = sums[parent] - sums[left];
        size_t lo = 0, hi = hash_count;
        while (lo < hi) {
            size_t mid = (lo + hi) / 2;
            if (hashes[mid].sum < target)
                lo = mid + 1;
            else
                hi = mid;
        }
        for (size_t j = lo; j < hash_count && hashes[j].sum == target; ++j) {
            Id right = hashes[j].id;
            if (hashes[j].count != counts[parent] - counts[left])
                continue;
            GSet c = sets[right];
            bool exact = true;
            uint32_t end = max(a.words, c.words);
            for (uint32_t w = 0; w < end; ++w) {
                U64 residual = (w < a.words ? words[a.offset + w] : 0ULL) &
                               ~(w < b.words ? words[b.offset + w] : 0ULL);
                if (residual != (w < c.words ? words[c.offset + w] : 0ULL)) {
                    exact = false;
                    break;
                }
            }
            if (exact) {
                U64 slot = atomicAdd(emitted, 1ULL);
                result[slot] = {parent, left, right};
                break;
            }
        }
    }
}
U64 cross_mix(U64 x) {
    x ^= x >> 30;
    x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27;
    x *= 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
} // namespace
std::vector<Split> cuda_cross_candidates(const CandidateSpace &space, const Config &cfg) {
    std::string reason;
    if (!cuda_available(&reason))
        throw std::runtime_error(reason);
    std::vector<GSet> sets;
    std::vector<U64> words(space.sets.arena.begin(), space.sets.arena.end());
    std::vector<uint32_t> counts;
    std::vector<U64> sums;
    std::vector<Id> parents, sides;
    std::vector<CrossHash> hashes;
    const U64 seed = cross_mix(0xDEADBEEFCAFEULL);
    for (Id i = 0; i < space.sets.records.size(); ++i) {
        const auto &r = space.sets.records[i];
        sets.push_back({r.offset, r.words});
        counts.push_back(r.count);
        U64 sum = 0;
        for (Id tax : members(space.sets.bits(i))) {
            U64 value = cross_mix(tax + seed);
            sum += value ? value : 1;
        }
        sums.push_back(sum);
        if (i < space.eligible.size() && space.eligible[i]) {
            parents.push_back(i);
            sides.push_back(i);
            hashes.push_back({sum, i, r.count});
        }
    }
    if (sides.empty())
        return {};
    parents.push_back(space.root);
    std::sort(sides.begin(), sides.end(),
              [&](Id a, Id b) { return counts[a] != counts[b] ? counts[a] < counts[b] : a < b; });
    std::sort(hashes.begin(), hashes.end(),
              [](auto &a, auto &b) { return a.sum != b.sum ? a.sum < b.sum : a.id < b.id; });
    size_t static_bytes = sets.size() * sizeof(GSet) + words.size() * sizeof(U64) +
                          counts.size() * sizeof(uint32_t) + sums.size() * sizeof(U64) +
                          parents.size() * sizeof(Id) + sides.size() * sizeof(Id) +
                          hashes.size() * sizeof(CrossHash);
    size_t available = budget(cfg);
    if (static_bytes > available * 3 / 4)
        throw std::runtime_error("cross-tree static species index exceeds CUDA memory budget");
    DeviceArray<GSet> ds;
    DeviceArray<U64> dw, dh_sums, counter;
    DeviceArray<uint32_t> dc;
    DeviceArray<Id> dp, db;
    DeviceArray<CrossHash> dh;
    ds.upload(sets);
    dw.upload(words);
    dc.upload(counts);
    dh_sums.upload(sums);
    dp.upload(parents);
    db.upload(sides);
    dh.upload(hashes);
    size_t capacity = std::max<size_t>(
        1, std::min<size_t>(cfg.gpu_dp_cap, (available - static_bytes) / sizeof(CrossResult) / 2));
    const size_t parent_batch = std::min<size_t>(64, capacity);
    const size_t side_batch = std::max<size_t>(1, std::min<size_t>(8192, capacity / parent_batch));
    DeviceArray<CrossResult> output(parent_batch * side_batch);
    counter.reserve(1);
    std::vector<CrossResult> host(parent_batch * side_batch);
    std::vector<Split> result;
    auto started = std::chrono::steady_clock::now();
    double last = 0;
    int progress_steps = 0;
    for (size_t pb = 0; pb < parents.size(); pb += parent_batch) {
        size_t pc = std::min(parent_batch, parents.size() - pb);
        for (size_t sb = 0; sb < sides.size(); sb += side_batch) {
            size_t se = std::min(sides.size(), sb + side_batch);
            counter.clear(1);
            cross_kernel<<<static_cast<unsigned>(pc), block_size>>>(
                ds.data(), dw.data(), dh.data(), hashes.size(), dp.data(), pb, pc, db.data(),
                dc.data(), dh_sums.data(), sb, se, output.data(), counter.data());
            check(cudaGetLastError(), "launch cross-tree candidates");
            U64 found = 0;
            counter.download(&found, 1);
            if (found > host.size())
                throw std::runtime_error("cross-tree append exceeded proven tile bound");
            output.download(host.data(), found);
            for (size_t j = 0; j < found; ++j)
                result.push_back({host[j].parent, host[j].left, host[j].right});
        }
        double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        if (cfg.gpu_dp_progress > 0 && elapsed - last >= cfg.gpu_dp_progress &&
            progress_steps < cfg.gpu_dp_steps) {
            log(cfg,
                "CUDA cross-tree parents: " + std::to_string(std::min(parents.size(), pb + pc)) +
                    "/" + std::to_string(parents.size()));
            last = elapsed;
            ++progress_steps;
        }
    }
    std::sort(result.begin(), result.end(), [&](const Split &a, const Split &b) {
        if (a.parent != b.parent) {
            if (a.parent == space.root)
                return false;
            if (b.parent == space.root)
                return true;
            return a.parent < b.parent;
        }
        if (counts[a.left] != counts[b.left])
            return counts[a.left] < counts[b.left];
        return a.left < b.left;
    });
    log(cfg, "CUDA cross-tree construction: " + std::to_string(result.size()) +
                 " candidate transitions");
    return result;
}
} // namespace stelar
