#include "stelar/core.hpp"
#include <iostream>
#include <numeric>
#include <random>
using namespace stelar;
static void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}

static uint64_t brute(const Bits &a, const Bits &b, const std::vector<Bits> &children) {
    uint64_t result = 0;
    for (size_t sister = 0; sister < children.size(); ++sister) {
        const auto x = members(children[sister]);
        for (size_t outside = 0; outside < children.size(); ++outside)
            if (outside != sister)
                for (Id z : members(children[outside]))
                    for (size_t i = 0; i < x.size(); ++i)
                        for (size_t j = i + 1; j < x.size(); ++j)
                            result += (contains(a, x[i]) && contains(a, x[j]) && contains(b, z)) ||
                                      (contains(b, x[i]) && contains(b, x[j]) && contains(a, z));
    }
    return result;
}

int main() {
    try {
        Registry quoted;
        auto parsed = parse_newick("[outer[nested]](('a b':.2,'it''s':.4)S,C)D;", quoted);
        check(quoted.size() == 3 && quoted.find("a b") != none && quoted.find("it's") != none,
              "quoted label or nested comment parsing");
        check(parsed.nodes[parsed.root].event == DUPLICATION,
              "root event tag must survive parsing");
        Registry roundtrip;
        check(parse_newick(newick(parsed, quoted, true), roundtrip).leaves.size() == 3,
              "quoted label round trip");
        bool rejected = false;
        try {
            parse_newick("((A,B),C;", quoted);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        check(rejected, "malformed Newick must be rejected");
        std::atomic<size_t> tasks{0};
        for (int threads : {4, 2, 1, 4})
            parallel_for(23, threads,
                         [&](size_t) { parallel_for(7, 4, [&](size_t) { ++tasks; }); });
        check(tasks == 4 * 23 * 7, "persistent pool reconfiguration/nested tasks");
        rejected = false;
        try {
            parallel_for(10, 3, [](size_t i) {
                if (i == 3)
                    throw std::invalid_argument("intentional worker exception");
            });
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        check(rejected, "worker exceptions must propagate without deadlocking");
        tasks = 0;
        parallel_for(19, 3, [&](size_t) { ++tasks; });
        check(tasks == 19, "pool remains usable after a failed task");
        SetPool pool;
        const Id a = pool.intern({1, 2, 0, 0});
        check(a == pool.intern({1, 2}), "set interning must trim trailing zero words");
        check(pool.count(a) == 2, "set cardinality");
        check(pool.intersect(a, pool, a) == 2, "set intersection");
        check(pool.find({1, 2, 0}) == a, "set lookup zero trimming");
        CandidateSpace space(130);
        const Id p = space.add_cluster(set_union(singleton(0), singleton(129)));
        const Id l = space.add_cluster(singleton(0)), r = space.add_cluster(singleton(129));
        space.add_split(p, l, r);
        space.add_split(p, r, l);
        space.add_split(p, l, r);
        check(space.splits.size() == 1, "unordered child split deduplication");

        std::mt19937 rng(19137);
        for (int trial = 0; trial < 80; ++trial) {
            const int n = trial % 3 == 0 ? 130 : 17, degree = trial % 4 + 2;
            Bits left, right;
            std::vector<Bits> children(degree);
            for (Id t = 0; t < Id(n); ++t) {
                if (rng() % 3 == 0)
                    left = set_union(left, singleton(t));
                else if (rng() % 2)
                    right = set_union(right, singleton(t));
                if (rng() % 5)
                    children[rng() % degree] = set_union(children[rng() % degree], Bits{});
            }
            // Distinct child species sets, including empty/unobserved species.
            for (Id t = 0; t < Id(n); ++t)
                if (rng() % 5) {
                    const auto index = rng() % degree;
                    children[index] = set_union(children[index], singleton(t));
                }
            SetPool candidates;
            const Id li = candidates.intern(left), ri = candidates.intern(right);
            Observations obs;
            obs.partitions.push_back({0, uint32_t(degree), 7});
            for (const auto &bits : children)
                obs.child_sets.push_back(obs.sets.intern(bits));
            check(split_weight(candidates, li, ri, obs) == Score(brute(left, right, children)) * 7,
                  "binary/polytomy split formula disagrees with triplet enumeration");
        }

        Registry registry;
        auto gene = parse_newick("(((A,B),(C,D)),((A,B),(C,D)))D;", registry, true, false);
        auto observations = build_observations({gene});
        check(observations.partitions.size() == 3,
              "identical observed child partitions must deduplicate");
        for (const auto &observation : observations.partitions)
            check(observation.frequency == 2, "deduplicated observation frequency");
        SetPool taxa;
        Id ab =
            taxa.intern(set_union(singleton(registry.find("A")), singleton(registry.find("B"))));
        Id cd =
            taxa.intern(set_union(singleton(registry.find("C")), singleton(registry.find("D"))));
        check(split_weight(taxa, ab, cd, observations) == 8,
              "multicopy root-event exact triplet score");

        // The active exact type must remain exact beyond unsigned 64-bit scores.
        Bits wide_a, wide_b;
        for (Id i = 0; i < 10000; ++i) {
            auto &bits = i < 5000 ? wide_a : wide_b;
            if (bits.size() <= i / 64)
                bits.resize(i / 64 + 1);
            bits[i / 64] |= uint64_t(1) << (i % 64);
        }
        SetPool wide;
        const Id wa = wide.intern(wide_a), wb = wide.intern(wide_b);
        Observations large;
        large.child_sets = {large.sets.intern(wide_a), large.sets.intern(wide_b)};
        large.partitions.push_back({0, 2, UINT32_MAX});
        const Score expected = Score(5000) * 5000 * 9998 / 2 * UINT32_MAX;
        check(expected > Score(UINT64_MAX), "128-bit test must exceed uint64 range");
        check(split_weight(wide, wa, wb, large) == expected, "128-bit exact score overflow");
        std::cout << "Core compact sets/deduplication, binary/polytomy brute triplets, multicopy "
                     "tags and uint128 PASS\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Core test failed: " << error.what() << '\n';
        return 1;
    }
}
