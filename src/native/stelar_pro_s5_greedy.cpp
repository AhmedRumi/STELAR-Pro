// S5 reuses ASTRAL-Pro's placement/NNI implementation in-process. The
// resulting trees are candidate proposals; STELAR-Pro never uses its quartet DP.
#include <jni.h>

// The bundled source is a single translation unit.  Rename its CLI entry
// point so that this JNI library can call its search classes directly.
#define main astralProCliMainUnusedByS5
#include "../../ASTER-Linux/src/astral-pro.cpp"
#undef main

#ifdef STELAR_S5_GPU_BUILD
std::unique_ptr<QuartetScoreBackend> makeS5GpuBackend(
    const TripartitionInitializer& genes, int taxa);
#endif

namespace {

void throwJava(JNIEnv* env, const char* type, const string& message) {
    jclass cls = env->FindClass(type);
    if (cls != nullptr) env->ThrowNew(cls, message.c_str());
}

// A JVM may run several tests in one process.  ASTER's CLI uses global input
// buffers, so clear them before each in-process search invocation.
void resetInputState(int threads, int partitions, unsigned long long seed) {
    TEXT.clear();
    pos = 0;
    nodecnt = 0;
    K = 0;
    duploss = 0;
    rootNtagTrees.clear();
    leafname_mapping.clear();
    name2id.clear();
    id2name.clear();
    tripInit.nodes.clear();
    tripInit.leafParent.clear();
    tripInit.nodes.resize(partitions);
    tripInit.leafParent.resize(partitions);
    meta.nThreads = threads;
    resolvePolytomies = true;
    rootNtag = false;
    LOG.enabled = false;
    srand(static_cast<unsigned int>(seed));
    RND_GENERATOR.seed(seed);
}

int countTrees(const char* path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open S5 gene-tree input");
    return static_cast<int>(std::count(std::istreambuf_iterator<char>(input),
                                       std::istreambuf_iterator<char>(), ';'));
}

}  // namespace

#ifdef STELAR_S5_GPU_BUILD
#define S5_JNI_ENTRY Java_stelarx_pro_AProGreedySearch_generateGPU
#else
#define S5_JNI_ENTRY Java_stelarx_pro_AProGreedySearch_generateCPU
#endif

extern "C" JNIEXPORT jobjectArray JNICALL
S5_JNI_ENTRY(
    JNIEnv* env, jclass, jstring inputPath, jint threads, jlong seed, jint rounds) {
    if (inputPath == nullptr || threads < 1 || rounds < 1) {
        throwJava(env, "java/lang/IllegalArgumentException",
                  "S5 requires an input path, positive thread count and positive round count");
        return nullptr;
    }
    const char* path = env->GetStringUTFChars(inputPath, nullptr);
    if (path == nullptr) return nullptr;
    try {
        const int geneCount = countTrees(path);
#ifdef STELAR_S5_GPU_BUILD
        // A GPU worker owns one gene family, so its color counter is race-free.
        const int partitions = geneCount;
#else
        // CPU workers own disjoint groups of gene families.
        const int partitions = threads;
#endif
        if (partitions < 1) throw std::invalid_argument("S5 input contains no trees");
        resetInputState(threads, partitions, static_cast<unsigned long long>(seed));
        annotate(path, "");
        env->ReleaseStringUTFChars(inputPath, path);
        path = nullptr;
        if (K == 0 || id2name.size() < 4) {
            throwJava(env, "java/lang/IllegalArgumentException",
                      "S5 requires at least one gene tree and four species");
            return nullptr;
        }

        ThreadPool pool;
        pool.initialize(threads);
        const int n = static_cast<int>(id2name.size());
        const int maxNni = 20 + static_cast<int>(2 * sqrt(n) * log2(n));
        ConstrainedOptimizationAlgorithm search(n, tripInit, id2name, pool, maxNni);
        vector<string> greedy;
        greedy.reserve(rounds);
        for (int round = 0; round < rounds; ++round) {
#ifdef STELAR_S5_GPU_BUILD
            auto backend = makeS5GpuBackend(tripInit, n);
            PlacementAlgorithm tree(search.taxonHash, tripInit, pool, maxNni,
                                    backend.get());
#else
            PlacementAlgorithm tree(search.taxonHash, tripInit, pool, maxNni);
#endif
            // roundId remains zero: independent random-order greedy searches,
            // without A-Pro's final quartet DP or its quartet-backed subsampling.
            search.createPlacementAlgorithm(tree, 0.0);
            tree.run();
            greedy.push_back(tree.printTree(id2name));
        }

        jclass stringClass = env->FindClass("java/lang/String");
        if (stringClass == nullptr) return nullptr;
        jobjectArray result = env->NewObjectArray(rounds, stringClass, nullptr);
        if (result == nullptr) return nullptr;
        for (int i = 0; i < rounds; ++i) {
            jstring newick = env->NewStringUTF(greedy[i].c_str());
            if (newick == nullptr) return nullptr;
            env->SetObjectArrayElement(result, i, newick);
            env->DeleteLocalRef(newick);
        }
        return result;
    } catch (const std::exception& e) {
        if (path != nullptr) env->ReleaseStringUTFChars(inputPath, path);
        throwJava(env, "java/lang/IllegalStateException", e.what());
        return nullptr;
    }
}
