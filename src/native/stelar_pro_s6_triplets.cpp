// S6 CPU scoring: a dense tensor remains for parity tests, while production
// uses compact speciation-child memberships and mutable colour states.
#include <jni.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>
#ifdef __GLIBC__
#include <malloc.h>
#endif

extern "C" JNIEXPORT jlongArray JNICALL
Java_stelarx_pro_S6TripletSupport_buildNative(JNIEnv* env, jclass, jint n,
    jintArray jPartOffsets, jintArray jChildOffsets, jintArray jTaxa, jint threads) {
    try {
        if (n < 3 || threads < 1) throw std::invalid_argument("invalid S6 dimensions");
        const jsize partCount = env->GetArrayLength(jPartOffsets) - 1;
        const jsize childCount = env->GetArrayLength(jChildOffsets) - 1;
        const jsize taxonCount = env->GetArrayLength(jTaxa);
        std::vector<jint> parts(partCount + 1), children(childCount + 1), taxa(taxonCount);
        env->GetIntArrayRegion(jPartOffsets, 0, partCount + 1, parts.data());
        env->GetIntArrayRegion(jChildOffsets, 0, childCount + 1, children.data());
        env->GetIntArrayRegion(jTaxa, 0, taxonCount, taxa.data());
        if (env->ExceptionCheck()) return nullptr;
        const int64_t cells = int64_t(n) * n * n;
        if (cells > 2147483647) throw std::invalid_argument("S6 tensor too large");
        std::vector<jlong> result(static_cast<size_t>(cells), 0);
        std::atomic<int> next{0};
        const int workers = std::min<int>(threads, n);
        std::vector<std::thread> pool;
        pool.reserve(workers);
        for (int worker = 0; worker < workers; worker++) {
            pool.emplace_back([&] {
                for (int c; (c = next.fetch_add(1)) < n;) {
                    jlong* row = result.data() + static_cast<size_t>(c) * n * n;
                    for (int p = 0; p < partCount; p++) {
                        for (int out = parts[p]; out < parts[p + 1]; out++) {
                            if (!std::binary_search(taxa.begin() + children[out],
                                    taxa.begin() + children[out + 1], c)) continue;
                            for (int side = parts[p]; side < parts[p + 1]; side++) {
                                if (side == out) continue;
                                const int end = children[side + 1];
                                for (int i = children[side]; i < end; i++) {
                                    const int a = taxa[i];
                                    if (a == c) continue;
                                    for (int j = i + 1; j < end; j++) {
                                        const int b = taxa[j];
                                        if (b != c) ++row[a * n + b];
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }
        for (auto& worker : pool) worker.join();
        jlongArray output = env->NewLongArray(static_cast<jsize>(cells));
        if (output == nullptr) return nullptr;
        env->SetLongArrayRegion(output, 0, static_cast<jsize>(cells), result.data());
        return output;
    } catch (const std::exception& error) {
        jclass cls = env->FindClass("java/lang/IllegalStateException");
        if (cls) env->ThrowNew(cls, error.what());
        return nullptr;
    }
}

namespace {

struct ColourIndex {
    struct State {
        std::vector<uint16_t> counts;
        std::vector<int8_t> colours;
        std::vector<int> desired, marks, touched;
        std::array<int64_t, 3> totals{0, 0, 0};
        int epoch = 0;

        State(int taxa, int children, int parts)
            : counts(static_cast<size_t>(3) * children, 0),
              colours(taxa, -1), desired(taxa, -1),
              marks(parts, 0) { touched.reserve(parts); }
    };

    int n;
    std::vector<int> parts, partOfChild, offsets, memberships;
    std::vector<State> states;

    ColourIndex(int taxa, int stateCount, std::vector<int> p,
                std::vector<int> childPart, std::vector<int> taxonOffsets,
                std::vector<int> taxonChildren)
        : n(taxa), parts(std::move(p)), partOfChild(std::move(childPart)),
          offsets(std::move(taxonOffsets)), memberships(std::move(taxonChildren)) {
        states.reserve(stateCount);
        for (int i = 0; i < stateCount; ++i)
            states.emplace_back(n, static_cast<int>(partOfChild.size()),
                                static_cast<int>(parts.size()) - 1);
    }

    void accumulate(State& state, int part, int sign) const {
        int64_t sumA = 0, sumB = 0, sumC = 0;
        for (int child = parts[part]; child < parts[part + 1]; ++child) {
            int base = 3 * child;
            sumA += state.counts[base];
            sumB += state.counts[base + 1];
            sumC += state.counts[base + 2];
        }
        for (int child = parts[part]; child < parts[part + 1]; ++child) {
            int base = 3 * child;
            int64_t a = state.counts[base], b = state.counts[base + 1];
            int64_t c = state.counts[base + 2];
            state.totals[0] += sign * a * b * (sumC - c);
            state.totals[1] += sign * a * c * (sumB - b);
            state.totals[2] += sign * b * c * (sumA - a);
        }
    }

    void update(State& state) const {
        if (++state.epoch == std::numeric_limits<int>::max()) {
            std::fill(state.marks.begin(), state.marks.end(), 0);
            state.epoch = 1;
        }
        state.touched.clear();
        for (int taxon = 0; taxon < n; ++taxon) {
            if (state.desired[taxon] == state.colours[taxon]) continue;
            for (int at = offsets[taxon]; at < offsets[taxon + 1]; ++at) {
                int part = partOfChild[memberships[at]];
                if (state.marks[part] != state.epoch) {
                    state.marks[part] = state.epoch;
                    state.touched.push_back(part);
                }
            }
        }
        for (int part : state.touched) accumulate(state, part, -1);
        for (int taxon = 0; taxon < n; ++taxon) {
            int before = state.colours[taxon], after = state.desired[taxon];
            if (before == after) continue;
            for (int at = offsets[taxon]; at < offsets[taxon + 1]; ++at) {
                int base = 3 * memberships[at];
                if (before >= 0) --state.counts[base + before];
                if (after >= 0) ++state.counts[base + after];
            }
            state.colours[taxon] = static_cast<int8_t>(after);
        }
        for (int part : state.touched) accumulate(state, part, 1);
    }

    int64_t splitWeight(State& state) const {
        update(state);
        int64_t score = 0;
        for (size_t p = 0; p + 1 < parts.size(); ++p) {
            int64_t sumA = 0, sumB = 0;
            for (int child = parts[p]; child < parts[p + 1]; ++child) {
                sumA += state.counts[3 * child];
                sumB += state.counts[3 * child + 1];
            }
            for (int child = parts[p]; child < parts[p + 1]; ++child) {
                int64_t a = state.counts[3 * child];
                int64_t b = state.counts[3 * child + 1];
                score += a * (a - 1) / 2 * (sumB - b)
                       + b * (b - 1) / 2 * (sumA - a);
            }
        }
        return score;
    }
};

void throwColourError(JNIEnv* env, const std::exception& error) {
    jclass cls = env->FindClass("java/lang/IllegalStateException");
    if (cls) env->ThrowNew(cls, error.what());
}

ColourIndex& getColourIndex(jlong handle, jint stateId) {
    auto* index = reinterpret_cast<ColourIndex*>(handle);
    if (!index || stateId < 0 || stateId >= static_cast<jint>(index->states.size()))
        throw std::invalid_argument("invalid S6 native colour state");
    return *index;
}

}  // namespace

extern "C" JNIEXPORT jlong JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_createColourIndex(JNIEnv* env, jclass,
    jint n, jintArray jParts, jintArray jPartOfChild,
    jintArray jOffsets, jintArray jMemberships, jint stateCount) {
    try {
        int pLen = env->GetArrayLength(jParts);
        int cLen = env->GetArrayLength(jPartOfChild);
        int mLen = env->GetArrayLength(jMemberships);
        if (n < 3 || n > 32767 || stateCount < 1 || pLen < 1 ||
                env->GetArrayLength(jOffsets) != n + 1)
            throw std::invalid_argument("invalid S6 native colour dimensions");
        std::vector<int> parts(pLen), childPart(cLen), offsets(n + 1), members(mLen);
        env->GetIntArrayRegion(jParts, 0, pLen, parts.data());
        env->GetIntArrayRegion(jPartOfChild, 0, cLen, childPart.data());
        env->GetIntArrayRegion(jOffsets, 0, n + 1, offsets.data());
        env->GetIntArrayRegion(jMemberships, 0, mLen, members.data());
        if (env->ExceptionCheck()) return 0;
        auto* index = new ColourIndex(n, stateCount, std::move(parts),
            std::move(childPart), std::move(offsets), std::move(members));
        return reinterpret_cast<jlong>(index);
    } catch (const std::exception& error) {
        throwColourError(env, error);
        return 0;
    }
}

extern "C" JNIEXPORT void JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_scoreGroupsNative(JNIEnv* env, jclass,
    jlong handle, jint stateId, jintArray jDesired, jlongArray jResult) {
    try {
        ColourIndex& index = getColourIndex(handle, stateId);
        if (env->GetArrayLength(jDesired) != index.n || env->GetArrayLength(jResult) != 3)
            throw std::invalid_argument("invalid S6 group-score buffers");
        auto& state = index.states[stateId];
        env->GetIntArrayRegion(jDesired, 0, index.n, state.desired.data());
        if (env->ExceptionCheck()) return;
        index.update(state);
        jlong result[3] = {state.totals[0], state.totals[1], state.totals[2]};
        env->SetLongArrayRegion(jResult, 0, 3, result);
    } catch (const std::exception& error) {
        throwColourError(env, error);
    }
}

extern "C" JNIEXPORT jlong JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_scoreSplitNative(JNIEnv* env, jclass,
    jlong handle, jint stateId, jintArray jDesired) {
    try {
        ColourIndex& index = getColourIndex(handle, stateId);
        if (env->GetArrayLength(jDesired) != index.n)
            throw std::invalid_argument("invalid S6 split-score buffer");
        auto& state = index.states[stateId];
        env->GetIntArrayRegion(jDesired, 0, index.n, state.desired.data());
        if (env->ExceptionCheck()) return 0;
        return index.splitWeight(state);
    } catch (const std::exception& error) {
        throwColourError(env, error);
        return 0;
    }
}

extern "C" JNIEXPORT void JNICALL
Java_stelarx_pro_S6DynamicTripletSupport_destroyColourIndex(JNIEnv*, jclass,
    jlong handle) {
    auto* index = reinterpret_cast<ColourIndex*>(handle);
    bool large = index != nullptr && index->partOfChild.size() > 100000;
    delete index;
#ifdef __GLIBC__
    // The one-pass search releases a large temporary native index before DP.
    // Return its freed pages so the later weight stage does not inherit them.
    if (large) malloc_trim(0);
#endif
}
