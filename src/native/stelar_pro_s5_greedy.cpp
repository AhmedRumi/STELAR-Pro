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

jobjectArray stringsToJava(JNIEnv* env, const vector<string>& values) {
    jclass stringClass = env->FindClass("java/lang/String");
    if (stringClass == nullptr) return nullptr;
    jobjectArray out = env->NewObjectArray(values.size(), stringClass, nullptr);
    env->DeleteLocalRef(stringClass);
    if (out == nullptr) return nullptr;
    for (size_t i = 0; i < values.size(); ++i) {
        jstring value = env->NewStringUTF(values[i].c_str());
        if (value == nullptr) return nullptr;
        env->SetObjectArrayElement(out, i, value);
        env->DeleteLocalRef(value);
    }
    return out;
}

// Native A-Pro still builds subset guide topologies. Their backbone is inferred
// in Java by STELAR's existing rooted-triplet DP, never by A-Pro's quartet DP.
string resolveTripletBackbone(JNIEnv* env, jstring taggedInput,
                              const vector<string>& guides,
                              const vector<int>& selected,
                              const vector<string>& names, bool gpu) {
    vector<string> subsetNames;
    subsetNames.reserve(selected.size());
    for (int id : selected) subsetNames.push_back(names.at(id));
    jobjectArray guideArray = stringsToJava(env, guides);
    jobjectArray nameArray = stringsToJava(env, subsetNames);
    if (guideArray == nullptr || nameArray == nullptr) {
        throw std::runtime_error("S5 could not allocate backbone guide arrays");
    }
    jclass bridge = env->FindClass("stelarx/pro/AProGreedySearch");
    if (bridge == nullptr) throw std::runtime_error("S5 backbone bridge class not found");
    jmethodID method = env->GetStaticMethodID(bridge, "inferTripletBackbone",
        "(Ljava/lang/String;[Ljava/lang/String;[Ljava/lang/String;Z)Ljava/lang/String;");
    if (method == nullptr) throw std::runtime_error("S5 backbone bridge method not found");
    jstring inferred = static_cast<jstring>(env->CallStaticObjectMethod(
        bridge, method, taggedInput, guideArray, nameArray,
        gpu ? JNI_TRUE : JNI_FALSE));
    env->DeleteLocalRef(guideArray);
    env->DeleteLocalRef(nameArray);
    env->DeleteLocalRef(bridge);
    if (env->ExceptionCheck() || inferred == nullptr) {
        throw std::runtime_error("S5 rooted-triplet backbone inference failed");
    }
    const char* chars = env->GetStringUTFChars(inferred, nullptr);
    if (chars == nullptr) throw std::runtime_error("S5 backbone Newick is unavailable");
    string backbone(chars);
    env->ReleaseStringUTFChars(inferred, chars);
    env->DeleteLocalRef(inferred);
    return backbone;
}

// Later rounds reuse the previous full STELAR-DP topology, induced on the
// newly selected taxa. Java prunes and anchors it without quartet-DP scoring.
string resolvePreviousBackbone(JNIEnv* env, jstring previousTree,
                               const vector<int>& selected,
                               const vector<string>& names) {
    vector<string> subsetNames;
    subsetNames.reserve(selected.size());
    for (int id : selected) subsetNames.push_back(names.at(id));
    jobjectArray allNames = stringsToJava(env, names);
    jobjectArray selectedNames = stringsToJava(env, subsetNames);
    if (allNames == nullptr || selectedNames == nullptr)
        throw std::runtime_error("S5 could not allocate prior-backbone taxa");
    jclass bridge = env->FindClass("stelarx/pro/AProGreedySearch");
    if (bridge == nullptr) throw std::runtime_error("S5 prior-backbone bridge not found");
    jmethodID method = env->GetStaticMethodID(bridge, "inducePreviousBackbone",
        "(Ljava/lang/String;[Ljava/lang/String;[Ljava/lang/String;)Ljava/lang/String;");
    if (method == nullptr) throw std::runtime_error("S5 prior-backbone method not found");
    jstring inferred = static_cast<jstring>(env->CallStaticObjectMethod(
        bridge, method, previousTree, allNames, selectedNames));
    env->DeleteLocalRef(allNames);
    env->DeleteLocalRef(selectedNames);
    env->DeleteLocalRef(bridge);
    if (env->ExceptionCheck() || inferred == nullptr)
        throw std::runtime_error("S5 prior STELAR-DP tree could not be restricted");
    const char* chars = env->GetStringUTFChars(inferred, nullptr);
    if (chars == nullptr) throw std::runtime_error("S5 prior backbone is unavailable");
    string backbone(chars);
    env->ReleaseStringUTFChars(inferred, chars);
    env->DeleteLocalRef(inferred);
    return backbone;
}

}  // namespace

#ifdef STELAR_S5_GPU_BUILD
#define S5_JNI_ENTRY Java_stelarx_pro_AProGreedySearch_generateGPU
#else
#define S5_JNI_ENTRY Java_stelarx_pro_AProGreedySearch_generateCPU
#endif

extern "C" JNIEXPORT jobjectArray JNICALL
S5_JNI_ENTRY(
    JNIEnv* env, jclass, jstring inputPath, jstring taggedPath, jstring previousTree,
    jint threads, jlong seed, jint rounds) {
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
        if (n >= 100 && taggedPath == nullptr) {
            throw std::invalid_argument(
                "S5 two-step search requires the rooted/tagged scoring input");
        }
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
            if (n >= 100) {
                // Each round independently follows A-Pro's two-step subset
                // placement. Later stages restrict the prior STELAR-DP tree.
                auto resolver = [&](const vector<string>& guides,
                                    const vector<int>& selected) {
                    return resolveTripletBackbone(env, taggedPath, guides,
                        selected, id2name,
#ifdef STELAR_S5_GPU_BUILD
                        true
#else
                        false
#endif
                    );
                };
                if (previousTree == nullptr) search.twoStepWorkflow(tree, resolver);
                else {
                    auto prior = [&](const vector<int>& selected) {
                        return resolvePreviousBackbone(env, previousTree,
                            selected, id2name);
                    };
                    search.twoStepWorkflow(tree, {}, prior);
                }
            } else if (previousTree != nullptr) {
                // A-Pro's small-taxon subsampling uses about one-quarter of the
                // previous optimal tree, then places the remaining taxa.
                vector<int> order;
                for (int i = 0; i < n; ++i) order.push_back(i);
                shuffle(order.begin(), order.end(), RND_GENERATOR);
                const int subsetSize = std::min(n, std::max(4, n / 4));
                vector<int> selected(order.begin(), order.begin() + subsetSize);
                string backbone = resolvePreviousBackbone(env, previousTree,
                    selected, id2name);
                unordered_map<string, int> subsetIds;
                unordered_set<int> added;
                for (int id : selected) subsetIds[id2name[id]] = id;
                int cursor = 0;
                search.guideSubtree(tree, backbone, subsetIds, cursor, added);
                if ((int) added.size() != subsetSize || tree.rootNodeId < 0)
                    throw std::runtime_error("S5 prior small-taxon backbone is invalid");
                for (int i = subsetSize; i < n; ++i) tree.order.push_back(order[i]);
                tree.run();
            } else {
                search.createPlacementAlgorithm(tree, 0.0);
                tree.run();
            }
            greedy.push_back(tree.printTree(id2name));
        }

        return stringsToJava(env, greedy);
    } catch (const std::exception& e) {
        if (path != nullptr) env->ReleaseStringUTFChars(inputPath, path);
        if (!env->ExceptionCheck())
            throwJava(env, "java/lang/IllegalStateException", e.what());
        return nullptr;
    }
}
