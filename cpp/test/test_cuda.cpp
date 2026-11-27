#include "stelar/core.hpp"
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>

using namespace stelar;
namespace {
void ensure(bool value, const std::string &message) {
    if (!value)
        throw std::runtime_error(message);
}
void compare_weights(CandidateSpace &space, const Observations &observations, Config cfg) {
    std::vector<Id> ids(space.splits.size());
    std::iota(ids.begin(), ids.end(), 0);
    const auto values = cuda_weights(space, observations, ids, cfg);
    ensure(values.size() == ids.size(), "CUDA result length");
    for (size_t i = 0; i < ids.size(); ++i) {
        const auto &s = space.splits[i];
        const Score expected = split_weight(space.sets, s.left, s.right, observations);
        ensure(values[i] == expected, "CUDA weight mismatch at split " + std::to_string(i) + ": " +
                                          score_string(values[i]) + " versus " +
                                          score_string(expected));
    }
}
void small_and_multicopy(Config cfg) {
    Registry registry;
    std::vector<Tree> genes;
    for (const std::string input : {"(((a,b),(c,d)),(e,(f,g)));", "((a,b),(a,(c,d)))D;",
                                    "((a,b,c),((d,e),f),g);", "(((a,b,c),(d,e)),(f,g));"})
        genes.push_back(parse_newick(input, registry, true, false));
    // Frequency deduplication and polytomies must both agree with CPU scoring.
    const Tree first = genes[0], second = genes[1];
    genes.push_back(first);
    genes.push_back(second);
    auto observations = build_observations(genes);
    CandidateSpace space(registry.size());
    const uint32_t states = 2187; // every 3-colouring of seven species
    for (uint32_t state = 0; state < states; ++state) {
        Bits a, b;
        uint32_t s = state;
        for (Id x = 0; x < registry.size(); ++x) {
            uint32_t colour = s % 3;
            s /= 3;
            if (colour == 1)
                a = set_union(a, singleton(x));
            else if (colour == 2)
                b = set_union(b, singleton(x));
        }
        if (a.empty() || b.empty())
            continue;
        space.add_split(space.add_cluster(set_union(a, b)), space.add_cluster(a),
                        space.add_cluster(b));
    }
    cfg.gpu_batch_size = 73;
    compare_weights(space, observations, cfg);
    cfg.gpu_batch_size = 0;
    cfg.gpu_batches = 11;
    compare_weights(space, observations, cfg);
    cfg.gpu_batches = 0;
    cfg.no_gpu_batch = true;
    compare_weights(space, observations, cfg);
    cfg.no_gpu_batch = false;
    cfg.score_type = "int128";
    cfg.gpu_batch_size = 91;
    compare_weights(space, observations, cfg);
    auto focus = make_cuda_focus(observations, registry.size(), cfg);
    ensure(bool(focus), "CUDA focus unavailable");
    for (Id x = 0; x < registry.size(); ++x) {
        std::vector<Id> selected;
        for (Id y = registry.size(); y-- > 0;)
            if (y != x)
                selected.push_back(y);
        auto actual = focus->focus(x, selected),
             expected = cpu_focus(observations, registry.size(), x, selected);
        ensure(actual.pair == expected.pair, "CUDA unordered global pair orientation differs");
        ensure(actual.sister == expected.sister, "CUDA sister focus differs");
        // Test cache reuse as well as non-sorted selected species.
        auto reused = focus->focus(x, selected);
        ensure(reused.pair == actual.pair && reused.sister == actual.sister,
               "CUDA focused cache reuse");
        if (selected.size() > 3) {
            selected.resize(3);
            actual = focus->focus(x, selected);
            expected = cpu_focus(observations, registry.size(), x, selected);
            ensure(actual.pair == expected.pair && actual.sister == expected.sister,
                   "CUDA selected-cache replacement");
        }
    }
    std::cout << "CUDA 64/128 weights, batch controls, multicopy/polytomy frequency and focused "
                 "support parity passed\n";
}
void wide_score(Config cfg) {
    constexpr size_t n = 200000;
    CandidateSpace space(n);
    Bits a = full_set(n / 2), b = set_difference(full_set(n), a);
    Id left = space.add_cluster(a), right = space.add_cluster(b);
    space.add_split(space.root, left, right);
    Observations observations;
    observations.child_sets = {observations.sets.intern(a), observations.sets.intern(b)};
    observations.partitions.push_back({0, 2, UINT32_MAX});
    Score expected = split_weight(space.sets, left, right, observations);
    ensure(expected > Score(UINT64_MAX), "synthetic overflow oracle did not exceed uint64");
    compare_weights(space, observations, cfg);
    std::cout << "CUDA true >64-bit score parity passed (" << score_string(expected) << ")\n";
}
void similarity_and_boundaries(Config cfg) {
    Registry registry;
    for (int i = 0; i < 201; ++i)
        registry.add("t" + std::to_string(i));
    std::vector<Tree> components;
    for (const std::string input : {"(((t0,t63),(t64,t127)),(t128,(t129,t200)));",
                                    "((t0,t64,t129),(t63,(t128,t200)));", "((t127,t128),t200);"})
        components.push_back(parse_newick(input, registry, false, false));
    // Binary Euler RMQ and unresolved-polytomy LCA matrices both have exact
    // integer numerators/denominators and therefore bit-identical doubles.
    const Matrix expected = triplet_matrix(components, registry.size(), cfg, false);
    Matrix actual;
    cfg.gpu_tile = 3;
    ensure(cuda_triplet_matrix(components, registry.size(), cfg, actual),
           "CUDA matrix unavailable");
    ensure(actual.values == expected.values,
           "CUDA matrix differs from CPU at species-word boundaries");
    std::vector<Tree> genes = components;
    genes.push_back(parse_newick("((t0,t63),(t0,(t64,t200)))D;", registry, false, false));
    auto observations = build_observations(genes);
    CandidateSpace candidates(registry.size());
    for (const auto &t : components)
        candidates.add_tree(t);
    compare_weights(candidates, observations, cfg);
    auto focus = make_cuda_focus(observations, registry.size(), cfg);
    std::vector<Id> selected = {200, 129, 128, 127, 64, 63, 0};
    for (Id x : {Id(1), Id(64), Id(200)}) {
        auto sel = selected;
        sel.erase(std::remove(sel.begin(), sel.end(), x), sel.end());
        auto a = focus->focus(x, sel), e = cpu_focus(observations, registry.size(), x, sel);
        ensure(a.pair == e.pair && a.sister == e.sister,
               "CUDA focus boundary or absent-species mismatch");
    }
    std::cout << "CUDA similarity, absent species and 64/128/200 word-boundary parity passed\n";
}
void cross_tree(Config cfg) {
    CandidateSpace base(7);
    for (unsigned mask = 1; mask < 128; ++mask) {
        Bits set{mask};
        base.add_cluster(set);
    }
    CandidateSpace cpu = base, gpu = base;
    Config host = cfg;
    host.compute = Config::CPU;
    host.threads = 4;
    cpu.add_cross(host);
    cfg.gpu_dp_cap = 32;
    gpu.add_cross(cfg);
    ensure(cpu.splits.size() == gpu.splits.size(), "GPU cross-tree split count");
    for (size_t i = 0; i < cpu.splits.size(); ++i) {
        const auto &a = cpu.splits[i];
        const auto &b = gpu.splits[i];
        ensure(a.parent == b.parent && a.left == b.left && a.right == b.right,
               "GPU cross-tree exact IDs/order");
    }
    std::cout << "CUDA cross-tree exact subset/residual, bounded tiles and deterministic order "
                 "parity passed\n";
}
void streamed_focus(Config cfg) {
    // Sparse high-ID sets force the focus index past the 1 MiB resident cap
    // without retaining gene-tree objects or allocating an n^3 tensor.
    Observations o;
    Id ab = o.sets.intern(Bits{3});
    for (Id i = 0; i < 200; ++i) {
        Id out = o.sets.intern(singleton(190000 + i));
        o.partitions.push_back({o.child_sets.size(), 2, 3});
        o.child_sets.push_back(ab);
        o.child_sets.push_back(out);
    }
    cfg.gpu_treewalk_cap_mb = 1;
    auto gpu = make_cuda_focus(o, 200000, cfg);
    for (Id x : {Id(190000), Id(190050), Id(190199)}) {
        std::vector<Id> selected = {1, 0, 190025, 190099};
        auto a = gpu->focus(x, selected), b = cpu_focus(o, 200000, x, selected);
        ensure(a.pair == b.pair && a.sister == b.sister, "streamed GPU focus differs from CPU");
    }
    std::cout << "CUDA bounded nonresident focused-support staging parity passed\n";
}
} // namespace
int main() {
    try {
        std::string reason;
        if (!cuda_available(&reason)) {
            std::cerr << "SKIP: " << reason << '\n';
            return 77;
        }
        Config cfg;
        cfg.compute = Config::GPU;
        cfg.gpu_strict = true;
        cfg.quiet = true;
        cfg.threads = 4;
        small_and_multicopy(cfg);
        wide_score(cfg);
        similarity_and_boundaries(cfg);
        cross_tree(cfg);
        streamed_focus(cfg);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "CUDA parity failure: " << e.what() << '\n';
        return 1;
    }
}
