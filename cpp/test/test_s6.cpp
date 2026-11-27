// Standalone white-box parity tests. Include the translation unit so placement
// deltas and incremental NNI counts can be tested without a public test API.
#include "../src/s6.cpp"
#include <cassert>
#include <iostream>
#include <set>

using namespace stelar;

namespace {
Score fixed_score(const GreedyTree &tree, const Observations &observations) {
    SetPool sets;
    std::vector<Id> clusters(tree.nodes.size(), none);
    Score total = 0;
    for (Id id : tree.postorder) {
        const GreedyNode &node = tree.nodes[id];
        if (node.leaf())
            clusters[id] = sets.intern(singleton(node.taxon));
        else {
            clusters[id] = sets.intern(
                set_union(sets.bits(clusters[node.left]), sets.bits(clusters[node.right])));
            total += split_weight(sets, clusters[node.left], clusters[node.right], observations);
        }
    }
    return total;
}

uint64_t direct_triplet(const Observations &observed, Id a, Id b, Id c) {
    uint64_t result = 0;
    for (const Observation &part : observed.partitions)
        for (uint32_t i = 0; i < part.children; ++i) {
            Bits paired = observed.sets.bits(observed.child_sets[part.offset + i]);
            if (!contains(paired, a) || !contains(paired, b))
                continue;
            for (uint32_t j = 0; j < part.children; ++j)
                if (i != j && contains(observed.sets.bits(observed.child_sets[part.offset + j]), c))
                    result += part.frequency;
        }
    return result;
}

void check_deltas_and_colours() {
    Registry registry;
    std::vector<Tree> genes;
    genes.push_back(parse_newick("(((a,b)S,(a,c)S)D,((d,e)S,f)S)S;", registry));
    genes.push_back(parse_newick("((a,(b,c)S)S,((d,e)S,f)S)S;", registry));
    genes.push_back(genes.back());
    genes.push_back(parse_newick("((a,b,c)S,(d,e,f)S)S;", registry));
    Observations observations = build_observations(genes);
    std::vector<Id> identity(registry.size());
    std::iota(identity.begin(), identity.end(), 0);
    GreedyTree tree(registry.size());
    Tree initial = parse_newick("((a,b),c);", registry, false);
    tree.import(initial, identity);
    Id x = registry.find("d");
    Focus focused = cpu_focus(observations, registry.size(), x, tree.leaves);
    JavaRandom rng(12345);
    locate(tree, x, focused, rng);
    Score before = fixed_score(tree, observations);
    std::vector<std::pair<Id, Score>> stack{{tree.root, 0}};
    while (!stack.empty()) {
        auto [id, prefix] = stack.back();
        stack.pop_back();
        const GreedyNode &node = tree.nodes[id];
        Score predicted = prefix + node.pair;
        GreedyTree copy = tree;
        Id parent = copy.nodes[id].parent;
        Id leaf = copy.leaf(x), joint = copy.joint(id, leaf);
        copy.nodes[joint].parent = parent;
        if (parent == none)
            copy.root = joint;
        else if (copy.nodes[parent].left == id)
            copy.nodes[parent].left = joint;
        else
            copy.nodes[parent].right = joint;
        copy.refresh();
        assert(fixed_score(copy, observations) - before == predicted);
        if (!node.leaf()) {
            stack.emplace_back(node.right, prefix + node.right_delta + tree.nodes[node.left].pair);
            stack.emplace_back(node.left, prefix + node.left_delta + tree.nodes[node.right].pair);
        }
    }
    GreedyTree full(registry.size());
    full.import(parse_newick("(((a,b),c),(d,(e,f)));", registry, false), identity);
    ColourIndex index(observations, registry.size());
    ColourState<uint16_t> narrow(index);
    ColourState<uint32_t> wide(index);
    for (int repeat = 0; repeat < 10; ++repeat)
        for (Id id : full.postorder) {
            const GreedyNode &p = full.nodes[id];
            if (p.leaf())
                continue;
            for (int side = 0; side < 2; ++side) {
                Id child = side ? p.right : p.left, c = side ? p.left : p.right;
                if (full.nodes[child].leaf())
                    continue;
                Id a = full.nodes[child].left, b = full.nodes[child].right;
                GroupScores expected{{0, 0, 0}};
                for (uint32_t ai = full.nodes[a].lo; ai < full.nodes[a].hi; ++ai)
                    for (uint32_t bi = full.nodes[b].lo; bi < full.nodes[b].hi; ++bi)
                        for (uint32_t ci = full.nodes[c].lo; ci < full.nodes[c].hi; ++ci) {
                            Id ta = full.leaves[ai], tb = full.leaves[bi], tc = full.leaves[ci];
                            expected[0] += direct_triplet(observations, ta, tb, tc);
                            expected[1] += direct_triplet(observations, ta, tc, tb);
                            expected[2] += direct_triplet(observations, tb, tc, ta);
                        }
                assert(narrow.groups(full, a, b, c) == expected);
                assert(wide.groups(full, a, b, c) == expected);
            }
        }
}

void check_group_matrices() {
    Registry registry;
    std::vector<Tree> components;
    components.push_back(parse_newick("((a,b),(c,(d,(e,f))));", registry));
    components.push_back(parse_newick("(((a,e),d),(c,(b,f)));", registry));
    components.push_back(components.front());
    components.push_back(parse_newick("(a,c);", registry, false));
    std::vector<std::vector<Id>> groups{{0, 2, 4}, {1, 3, 5}};
    Config cfg;
    cfg.quiet = true;
    cfg.compute = Config::CPU;
    auto matrices = partition_matrices(components, groups, registry.size(), cfg);
    for (size_t g = 0; g < groups.size(); ++g) {
        std::vector<int> map(registry.size(), -1);
        for (Id i = 0; i < groups[g].size(); ++i)
            map[groups[g][i]] = i;
        std::vector<Tree> restricted;
        for (const Tree &source : components) {
            Tree tree = restrict_tree(source, map);
            if (tree.leaves.size() >= 3)
                restricted.push_back(std::move(tree));
        }
        Matrix reference = triplet_matrix(restricted, groups[g].size(), cfg, false);
        for (size_t i = 0; i < reference.values.size(); ++i)
            assert(reference.values[i] == matrices[g].values[i]);
    }
}

std::string balanced(Id first, Id end) {
    if (first + 1 == end)
        return "s" + std::to_string(first);
    Id middle = first + (end - first) / 2;
    return "(" + balanced(first, middle) + "," + balanced(middle, end) + ")";
}

void check_full_path() {
    Registry registry;
    std::vector<Tree> genes;
    genes.push_back(parse_newick(balanced(0, 129) + ";", registry));
    genes.push_back(genes.front());
    const std::string original = newick(genes.front(), registry, true);
    Config serial;
    serial.compute = Config::CPU;
    serial.quiet = true;
    serial.seed = 0xDEADBEEFCAFE;
    serial.rounds = 2;
    serial.threads = 1;
    auto one = generate_s6(genes, registry, serial);
    Config parallel = serial;
    parallel.threads = 3;
    auto three = generate_s6(genes, registry, parallel);
    assert(one.size() == 2 && three.size() == 2);
    for (size_t round = 0; round < one.size(); ++round) {
        assert(one[round].binary());
        assert(one[round].leaves.size() == registry.size());
        assert(cardinality(one[round].taxa()) == registry.size());
        assert(newick(one[round], registry) == newick(three[round], registry));
    }
    assert(newick(genes.front(), registry, true) == original);
    Config cross = serial;
    cross.rounds = 1;
    cross.cross_tree = true;
    auto recombined = generate_s6(genes, registry, cross);
    assert(recombined.size() == 1 && recombined[0].leaves.size() == registry.size());
}
} // namespace

int main() {
    check_deltas_and_colours();
    check_group_matrices();
    check_full_path();
    std::cout << "S3 tests passed: insertion deltas, multicopy/polytomy colour scores, "
                 "one-pass group matrix parity, >100-species rounds and CPU thread determinism\n";
}
