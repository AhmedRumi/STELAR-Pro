// Build the S6 rooted-triplet support tensor from original tagged gene nodes.
// Outgroup-major ownership gives every worker disjoint output, without atomics.
#include <jni.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <thread>
#include <vector>

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
