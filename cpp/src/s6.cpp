#include "stelar/core.hpp"

#include <array>
#include <climits>
#include <cmath>
#include <numeric>
#include <queue>
#include <sstream>
#include <type_traits>

namespace stelar {
namespace {

// Stable integer node identifiers replace the Java greedy-node graph. Each
// subtree refers to one interval in a shared DFS leaf ordering, rather than
// owning its own array of descendant species.
struct GreedyNode {
    Id taxon = none, left = none, right = none, parent = none;
    uint32_t lo = 0, hi = 0;
    Score pair = 0, left_delta = 0, right_delta = 0;
    bool leaf() const {
        return taxon != none;
    }
};

struct GreedyTree {
    std::vector<GreedyNode> nodes;
    std::vector<Id> leaves, postorder, preorder;
    Id root = none;

    explicit GreedyTree(size_t n) {
        nodes.reserve(2 * n);
        leaves.reserve(n);
    }

    Id leaf(Id taxon) {
        GreedyNode node;
        node.taxon = taxon;
        nodes.push_back(node);
        return static_cast<Id>(nodes.size() - 1);
    }

    Id joint(Id left, Id right) {
        GreedyNode node;
        node.left = left;
        node.right = right;
        nodes.push_back(node);
        Id id = static_cast<Id>(nodes.size() - 1);
        nodes[left].parent = nodes[right].parent = id;
        return id;
    }

    void refresh() {
        leaves.clear();
        postorder.clear();
        preorder.clear();
        std::vector<std::pair<Id, bool>> stack{{root, false}};
        while (!stack.empty()) {
            auto entry = stack.back();
            stack.pop_back();
            GreedyNode &node = nodes[entry.first];
            if (entry.second) {
                node.hi = static_cast<uint32_t>(leaves.size());
                postorder.push_back(entry.first);
                continue;
            }
            node.lo = static_cast<uint32_t>(leaves.size());
            preorder.push_back(entry.first);
            if (node.leaf()) {
                leaves.push_back(node.taxon);
                node.hi = static_cast<uint32_t>(leaves.size());
                postorder.push_back(entry.first);
            } else {
                stack.emplace_back(entry.first, true);
                stack.emplace_back(node.right, false);
                stack.emplace_back(node.left, false);
            }
        }
    }

    void import(const Tree &source, const std::vector<Id> &global_ids) {
        std::vector<Id> mapped(source.nodes.size(), none);
        for (Id id : source.postorder()) {
            const Node &node = source.nodes[id];
            if (node.child == none) {
                if (node.taxon >= global_ids.size())
                    throw std::runtime_error("Invalid S3 backbone species identifier");
                mapped[id] = leaf(global_ids[node.taxon]);
            } else {
                Id first = node.child;
                Id second = source.nodes[first].next;
                if (second == none || source.nodes[second].next != none)
                    throw std::runtime_error("S3 needs binary backbone candidates");
                mapped[id] = joint(mapped[first], mapped[second]);
            }
        }
        root = mapped[source.root];
        refresh();
    }

    Tree export_tree() const {
        Tree tree;
        tree.nodes.reserve(nodes.size());
        std::vector<Id> mapped(nodes.size(), none);
        for (Id id : postorder) {
            const GreedyNode &node = nodes[id];
            mapped[id] = node.leaf() ? tree.add_leaf(node.taxon)
                                     : tree.add_internal({mapped[node.left], mapped[node.right]},
                                                         SPECIATION);
        }
        tree.root = mapped[root];
        tree.index();
        return tree;
    }
};

Id locate(GreedyTree &tree, Id x, const Focus &focus, JavaRandom &random) {
    const size_t m = focus.selected.size();
    std::vector<Id> local(focus.n, none);
    for (Id i = 0; i < m; ++i)
        local[focus.selected[i]] = i;
    auto support = [&](Id a, Id b, Id c) -> uint64_t {
        if (a == b || a == c || b == c)
            return 0;
        if (c == x) {
            if (a > b)
                std::swap(a, b);
            return focus.pair[size_t(local[a]) * m + local[b]];
        }
        Id sister = a == x ? b : a;
        return focus.sister[size_t(local[c]) * m + local[sister]];
    };
    for (Id id : tree.postorder) {
        GreedyNode &node = tree.nodes[id];
        if (node.leaf()) {
            node.pair = 0;
            continue;
        }
        const GreedyNode &left = tree.nodes[node.left];
        const GreedyNode &right = tree.nodes[node.right];
        Score across = 0, left_delta = 0, right_delta = 0;
        for (uint32_t i = left.lo; i < left.hi; ++i) {
            Id a = tree.leaves[i];
            for (uint32_t j = right.lo; j < right.hi; ++j) {
                Id b = tree.leaves[j];
                across += support(a, b, x);
                left_delta += support(x, a, b);
                right_delta += support(x, b, a);
            }
        }
        node.pair = left.pair + right.pair + across;
        node.left_delta = left_delta;
        node.right_delta = right_delta;
    }
    Score best = 0;
    bool found = false;
    std::vector<Id> ties;
    std::vector<std::pair<Id, Score>> stack{{tree.root, 0}};
    while (!stack.empty()) {
        auto entry = stack.back();
        stack.pop_back();
        const GreedyNode &node = tree.nodes[entry.first];
        Score delta = entry.second + node.pair;
        if (!found || delta > best) {
            found = true;
            best = delta;
            ties.clear();
            ties.push_back(entry.first);
        } else if (delta == best)
            ties.push_back(entry.first);
        if (!node.leaf()) {
            stack.emplace_back(node.right,
                               entry.second + node.right_delta + tree.nodes[node.left].pair);
            stack.emplace_back(node.left,
                               entry.second + node.left_delta + tree.nodes[node.right].pair);
        }
    }
    if (ties.empty())
        throw std::runtime_error("S3 has no placement edge");
    return ties[random.next_int(static_cast<int>(ties.size()))];
}

// Immutable CSR child memberships are shared by the bounded mutable NNI states.
// Observation frequencies preserve all original-family multiplicities after
// partition deduplication; no triplet tensor or split-by-partition table exists.
struct ColourIndex {
    size_t taxa;
    std::vector<Id> parts, part_of_child, memberships;
    std::vector<uint64_t> offsets;
    std::vector<uint32_t> frequency;

    ColourIndex(const Observations &observed, size_t n) : taxa(n), offsets(n + 1, 0) {
        if (observed.child_sets.size() >= none)
            throw std::runtime_error("S3 child membership index exceeds 32-bit identifiers");
        parts.reserve(observed.partitions.size() + 1);
        frequency.reserve(observed.partitions.size());
        parts.push_back(0);
        part_of_child.resize(observed.child_sets.size());
        for (Id p = 0; p < observed.partitions.size(); ++p) {
            const Observation &part = observed.partitions[p];
            if (part.offset != parts.back())
                throw std::runtime_error("S3 observations must use packed child identifiers");
            for (uint32_t c = 0; c < part.children; ++c) {
                Id child = static_cast<Id>(part.offset + c);
                part_of_child[child] = p;
                Id set = observed.child_sets[child];
                const uint64_t *words = observed.sets.data(set);
                for (uint32_t w = 0; w < observed.sets.records[set].words; ++w) {
                    uint64_t word = words[w];
                    while (word) {
                        Id taxon = 64 * w + static_cast<Id>(__builtin_ctzll(word));
                        word &= word - 1;
                        if (taxon >= n)
                            throw std::runtime_error("Invalid S3 observation species");
                        ++offsets[taxon + 1];
                    }
                }
            }
            parts.push_back(static_cast<Id>(part.offset + part.children));
            frequency.push_back(part.frequency);
        }
        for (size_t i = 1; i < offsets.size(); ++i)
            offsets[i] += offsets[i - 1];
        memberships.resize(offsets.back());
        std::vector<uint64_t> cursor = offsets;
        for (Id child = 0; child < observed.child_sets.size(); ++child) {
            Id set = observed.child_sets[child];
            const uint64_t *words = observed.sets.data(set);
            for (uint32_t w = 0; w < observed.sets.records[set].words; ++w) {
                uint64_t word = words[w];
                while (word) {
                    Id taxon = 64 * w + static_cast<Id>(__builtin_ctzll(word));
                    word &= word - 1;
                    memberships[cursor[taxon]++] = child;
                }
            }
        }
    }
};

using SignedScore = __int128;
using GroupScores = std::array<SignedScore, 3>;

struct ColourStateBase {
    virtual ~ColourStateBase() = default;
    virtual GroupScores groups(const GreedyTree &, Id, Id, Id) = 0;
};

template <class Count> class ColourState final : public ColourStateBase {
    const ColourIndex &index;
    std::vector<Count> counts;
    std::vector<int8_t> colours, desired;
    std::vector<uint32_t> marks;
    std::vector<Id> touched;
    GroupScores totals{{0, 0, 0}};
    uint32_t epoch = 0;

    void accumulate(Id part, int sign) {
        SignedScore sa = 0, sb = 0, sc = 0;
        for (Id child = index.parts[part]; child < index.parts[part + 1]; ++child) {
            sa += counts[3ULL * child];
            sb += counts[3ULL * child + 1];
            sc += counts[3ULL * child + 2];
        }
        SignedScore factor = SignedScore(sign) * index.frequency[part];
        for (Id child = index.parts[part]; child < index.parts[part + 1]; ++child) {
            SignedScore a = counts[3ULL * child], b = counts[3ULL * child + 1];
            SignedScore c = counts[3ULL * child + 2];
            totals[0] += factor * a * b * (sc - c);
            totals[1] += factor * a * c * (sb - b);
            totals[2] += factor * b * c * (sa - a);
        }
    }

    void update() {
        if (++epoch == 0) {
            std::fill(marks.begin(), marks.end(), 0);
            epoch = 1;
        }
        touched.clear();
        for (Id taxon = 0; taxon < index.taxa; ++taxon) {
            if (desired[taxon] == colours[taxon])
                continue;
            for (uint64_t at = index.offsets[taxon]; at < index.offsets[taxon + 1]; ++at) {
                Id part = index.part_of_child[index.memberships[at]];
                if (marks[part] != epoch) {
                    marks[part] = epoch;
                    touched.push_back(part);
                }
            }
        }
        for (Id part : touched)
            accumulate(part, -1);
        for (Id taxon = 0; taxon < index.taxa; ++taxon) {
            int before = colours[taxon], after = desired[taxon];
            if (before == after)
                continue;
            for (uint64_t at = index.offsets[taxon]; at < index.offsets[taxon + 1]; ++at) {
                size_t base = 3ULL * index.memberships[at];
                if (before >= 0)
                    --counts[base + before];
                if (after >= 0)
                    ++counts[base + after];
            }
            colours[taxon] = static_cast<int8_t>(after);
        }
        for (Id part : touched)
            accumulate(part, 1);
    }

  public:
    explicit ColourState(const ColourIndex &ix)
        : index(ix), counts(3ULL * ix.part_of_child.size(), 0), colours(ix.taxa, -1),
          desired(ix.taxa, -1), marks(ix.frequency.size(), 0) {
        touched.reserve(ix.frequency.size());
    }

    GroupScores groups(const GreedyTree &tree, Id a, Id b, Id c) override {
        std::fill(desired.begin(), desired.end(), int8_t(-1));
        const std::array<Id, 3> groups{{a, b, c}};
        for (int colour = 0; colour < 3; ++colour) {
            const GreedyNode &node = tree.nodes[groups[colour]];
            for (uint32_t pos = node.lo; pos < node.hi; ++pos)
                desired[tree.leaves[pos]] = static_cast<int8_t>(colour);
        }
        update();
        return totals;
    }
};

int nni_choice(const GreedyTree &tree, Id parent, ColourStateBase &scorer) {
    const GreedyNode &node = tree.nodes[parent];
    if (node.leaf())
        return 0;
    for (int side = 0; side < 2; ++side) {
        Id child = side == 0 ? node.left : node.right;
        if (tree.nodes[child].leaf())
            continue;
        Id sibling = side == 0 ? node.right : node.left;
        GroupScores scores =
            scorer.groups(tree, tree.nodes[child].left, tree.nodes[child].right, sibling);
        if (scores[1] > scores[0] || scores[2] > scores[0])
            return 1 + 2 * side + (scores[1] >= scores[2] ? 0 : 1);
    }
    return 0;
}

void refine_nni(GreedyTree &tree, const ColourIndex &index, const Config &cfg) {
    const size_t n = index.taxa;
    const int limit =
        20 + static_cast<int>(2 * std::sqrt(double(n)) * std::log(double(n)) / std::log(2.0));
    const int workers = n >= 100 ? std::max(1, std::min(3, cfg.threads)) : 1;
    std::vector<std::unique_ptr<ColourStateBase>> scorers;
    for (int i = 0; i < workers; ++i) {
        if (n <= UINT16_MAX)
            scorers.emplace_back(new ColourState<uint16_t>(index));
        else
            scorers.emplace_back(new ColourState<uint32_t>(index));
    }
    std::vector<int> cache(tree.nodes.size(), INT_MIN);
    int applied = 0;
    for (int move = 0; move < limit; ++move) {
        tree.refresh();
        Id parent = none;
        int choice = 0;
        if (workers == 1) {
            for (Id node : tree.postorder) {
                if (cache[node] == INT_MIN)
                    cache[node] = nni_choice(tree, node, *scorers[0]);
                if (cache[node]) {
                    parent = node;
                    choice = cache[node];
                    break;
                }
            }
        } else {
            const size_t window = 24 * workers;
            for (size_t begin = 0; begin < tree.postorder.size(); begin += window) {
                const size_t end = std::min(begin + window, tree.postorder.size());
                const size_t chunk = (end - begin + workers - 1) / workers;
                std::vector<int> choices(end - begin);
                size_t known = end;
                for (size_t i = begin; i < end; ++i) {
                    int value = choices[i - begin] = cache[tree.postorder[i]];
                    if (value != INT_MIN && value != 0 && i < known)
                        known = i;
                }
                bool missing = false;
                for (size_t i = begin; i < known; ++i)
                    if (choices[i - begin] == INT_MIN) {
                        missing = true;
                        break;
                    }
                if (!missing) {
                    if (known < end) {
                        parent = tree.postorder[known];
                        choice = choices[known - begin];
                        break;
                    }
                    continue;
                }
                std::atomic<size_t> earliest{known};
                parallel_for(workers, workers, [&](size_t worker) {
                    size_t first = begin + worker * chunk;
                    size_t last = std::min(end, first + chunk);
                    for (size_t i = first; i < last && i < earliest.load(); ++i) {
                        if (choices[i - begin] != INT_MIN)
                            continue;
                        int candidate = nni_choice(tree, tree.postorder[i], *scorers[worker]);
                        choices[i - begin] = candidate;
                        if (candidate) {
                            size_t current = earliest.load();
                            while (i < current && !earliest.compare_exchange_weak(current, i)) {
                            }
                        }
                    }
                });
                for (size_t i = begin; i < end; ++i)
                    if (choices[i - begin] != INT_MIN)
                        cache[tree.postorder[i]] = choices[i - begin];
                size_t first = earliest.load();
                if (first < end) {
                    parent = tree.postorder[first];
                    choice = choices[first - begin];
                    break;
                }
            }
        }
        if (parent == none)
            break;
        GreedyNode &p = tree.nodes[parent];
        bool first = choice <= 2;
        Id child = first ? p.left : p.right;
        Id sibling = first ? p.right : p.left;
        Id a = tree.nodes[child].left, b = tree.nodes[child].right;
        Id grandparent = p.parent;
        Id moved = (choice & 1) ? b : a, kept = (choice & 1) ? a : b;
        tree.nodes[child].left = kept;
        tree.nodes[child].right = sibling;
        tree.nodes[kept].parent = tree.nodes[sibling].parent = child;
        if (first)
            p.right = moved;
        else
            p.left = moved;
        tree.nodes[moved].parent = parent;
        for (Id changed : {grandparent, parent, child, sibling, a, b})
            if (changed != none)
                cache[changed] = INT_MIN;
        ++applied;
    }
    tree.refresh();
    log(cfg, "S3 triplet NNI: " + std::to_string(applied) + " improving moves");
}

struct GroupPair {
    double score;
    Id a, b;
};
struct GroupPairLess {
    bool operator()(const GroupPair &x, const GroupPair &y) const {
        if (x.score != y.score)
            return x.score < y.score;
        if (x.a != y.a)
            return x.a < y.a;
        return x.b < y.b;
    }
};

Id group_upgma(GreedyTree &tree, const std::vector<Id> &group, const Matrix &matrix) {
    size_t m = group.size();
    if (m == 1)
        return tree.leaf(group[0]);
    if (m == 2) {
        Id left = tree.leaf(group[0]), right = tree.leaf(group[1]);
        return tree.joint(left, right);
    }
    size_t capacity = 2 * m - 1;
    Matrix similarity(capacity);
    std::vector<Id> roots(capacity), sizes(capacity, 0);
    std::vector<uint8_t> active(capacity, 0);
    std::priority_queue<GroupPair, std::vector<GroupPair>, GroupPairLess> nearest;
    for (Id a = 0; a < m; ++a) {
        roots[a] = tree.leaf(group[a]);
        sizes[a] = 1;
        active[a] = 1;
        for (Id b = 0; b < a; ++b) {
            double score = matrix.get(a, b);
            similarity.set(a, b, score);
            nearest.push({score, b, a});
        }
    }
    for (Id next = static_cast<Id>(m); next < capacity; ++next) {
        while (!nearest.empty() && (!active[nearest.top().a] || !active[nearest.top().b]))
            nearest.pop();
        if (nearest.empty())
            throw std::runtime_error("S3 group UPGMA exhausted pairs");
        GroupPair best = nearest.top();
        nearest.pop();
        Id a = best.a, b = best.b;
        roots[next] = tree.joint(roots[a], roots[b]);
        sizes[next] = sizes[a] + sizes[b];
        active[a] = active[b] = 0;
        for (Id k = 0; k < next; ++k) {
            if (!active[k])
                continue;
            double score =
                (sizes[a] * similarity.get(a, k) + sizes[b] * similarity.get(b, k)) / sizes[next];
            similarity.set(next, k, score);
            nearest.push({score, k, next});
        }
        active[next] = 1;
    }
    return roots.back();
}

// All groups are disjoint. Accumulate their matrices directly during one
// postorder traversal of each observed component. Small-to-large buckets hold
// only live group leaves; restricted Tree copies are never materialized here.
// At a pair's LCA, k - |descendants in its group| is precisely the triplet
// numerator of the induced (uncompleted) component.
using GroupBuckets = std::unordered_map<Id, std::vector<Id>>;

GroupBuckets accumulate_groups(const Tree &source, Id node, const std::vector<int> &group_of,
                               const std::vector<Id> &local_of, const std::vector<uint32_t> &total,
                               std::vector<Matrix> &numerator, std::vector<Matrix> &denominator) {
    const Node &original = source.nodes[node];
    if (original.child == none) {
        int group = group_of[original.taxon];
        if (group < 0 || total[group] < 3)
            return {};
        GroupBuckets one;
        one.emplace(static_cast<Id>(group), std::vector<Id>{local_of[original.taxon]});
        return one;
    }
    Id first = original.child, second = source.nodes[first].next;
    if (second == none || source.nodes[second].next != none)
        throw std::runtime_error("S3 subset-local components must be binary");
    GroupBuckets left =
        accumulate_groups(source, first, group_of, local_of, total, numerator, denominator);
    GroupBuckets right =
        accumulate_groups(source, second, group_of, local_of, total, numerator, denominator);
    if (left.size() < right.size())
        left.swap(right);
    for (auto &entry : right) {
        Id group = entry.first;
        auto found = left.find(group);
        if (found == left.end()) {
            left.emplace(group, std::move(entry.second));
            continue;
        }
        std::vector<Id> &a = found->second;
        std::vector<Id> &b = entry.second;
        double num = double(total[group] - a.size() - b.size());
        double den = double(total[group] - 2);
        Matrix &nums = numerator[group];
        Matrix &dens = denominator[group];
        for (Id i : a)
            for (Id j : b) {
                nums.set(i, j, nums.get(i, j) + num);
                dens.set(i, j, dens.get(i, j) + den);
            }
        if (a.size() < b.size())
            a.swap(b);
        a.insert(a.end(), b.begin(), b.end());
    }
    return left;
}

std::vector<Matrix> partition_matrices(const std::vector<Tree> &components,
                                       const std::vector<std::vector<Id>> &groups, size_t n,
                                       const Config &cfg) {
    std::vector<int> group_of(n, -1);
    std::vector<Id> local_of(n, none);
    std::vector<Matrix> matrices, denominators;
    matrices.reserve(groups.size());
    denominators.reserve(groups.size());
    for (Id g = 0; g < groups.size(); ++g) {
        size_t m = groups[g].size();
        matrices.emplace_back(m >= 3 ? m : 0);
        denominators.emplace_back(m >= 3 ? m : 0);
        if (m < 3)
            continue;
        for (Id i = 0; i < m; ++i) {
            Id taxon = groups[g][i];
            if (group_of[taxon] >= 0)
                throw std::runtime_error("Overlapping S3 partition groups");
            group_of[taxon] = static_cast<int>(g);
            local_of[taxon] = i;
        }
    }
    std::vector<uint32_t> total(groups.size(), 0);
    for (const Tree &source : components) {
        std::fill(total.begin(), total.end(), 0);
        for (Id taxon : source.leaves)
            if (group_of[taxon] >= 0)
                ++total[group_of[taxon]];
        accumulate_groups(source, source.root, group_of, local_of, total, matrices, denominators);
    }
    for (size_t g = 0; g < matrices.size(); ++g) {
        Matrix &num = matrices[g];
        const Matrix &den = denominators[g];
        for (size_t i = 0; i < num.values.size(); ++i)
            num.values[i] = den.values[i] > 0 ? num.values[i] / den.values[i] : 0;
    }
    log(cfg, "S3 partition matrices: " + std::to_string(groups.size()) +
                 " disjoint groups in one DISCO traversal");
    return matrices;
}

std::vector<Tree> backbone_guides(const std::vector<Tree> &components,
                                  const std::vector<Id> &selected, size_t global_taxa,
                                  const Config &cfg) {
    const size_t n = selected.size();
    std::vector<int> map(global_taxa, -1);
    for (Id i = 0; i < n; ++i)
        map[selected[i]] = static_cast<int>(i);
    std::vector<Tree> restricted;
    restricted.reserve(components.size());
    for (const Tree &source : components) {
        Tree tree = restrict_tree(source, map);
        if (tree.root == none)
            continue;
        tree.scoped = true;
        tree.missing.clear();
        for (Id local = 0; local < n; ++local)
            if (contains(source.missing, selected[local])) {
                if (tree.missing.size() <= local / 64)
                    tree.missing.resize(local / 64 + 1, 0);
                tree.missing[local / 64] |= uint64_t(1) << (local % 64);
            }
        if (tree.leaves.size() + cardinality(tree.missing) >= 3)
            restricted.push_back(std::move(tree));
    }
    Matrix matrix = triplet_matrix(restricted, n, cfg, n >= 64);
    bool missing = false;
    for (const Tree &tree : restricted)
        if (cardinality(tree.missing)) {
            missing = true;
            break;
        }
    std::vector<Id> rows = missing ? sorted_rows(matrix, cfg) : std::vector<Id>{};
    parallel_for(restricted.size(), cfg.threads, [&](size_t i) {
        if (cardinality(restricted[i].missing))
            complete_tree(restricted[i], rows, n, true);
    });
    std::vector<Tree> guides;
    guides.reserve(restricted.size() + 1);
    for (Tree &tree : restricted)
        if (tree.leaves.size() >= 3)
            guides.push_back(std::move(tree));
    guides.push_back(upgma(std::move(matrix), cfg));
    log(cfg, "S3 backbone: " + std::to_string(n) + " taxa, " + std::to_string(guides.size() - 1) +
                 " completed duplication-scope components + UPGMA");
    return guides;
}

} // namespace

std::vector<Tree> generate_s6(const std::vector<Tree> &genes, const Registry &registry,
                              const Config &cfg) {
    const size_t n = registry.size();
    if (n <= 100) {
        // The executable routes this case through S2's complete DP pipeline.
        // Keep the helper safe for callers that construct proposals directly.
        throw std::invalid_argument("S3 with at most 100 taxa must use the S2 fallback pipeline");
    }
    if (cfg.rounds < 1)
        throw std::invalid_argument("S3 needs at least one round");
    DiscoResult disco = decompose(genes, n, Scope::DUPLICATION, cfg);
    if (disco.trees.empty())
        throw std::runtime_error("S3 needs at least one DISCO component");
    Observations observed = build_observations(genes);
    ColourIndex colour(observed, n);
    std::unique_ptr<CudaFocus> gpu;
    if (cfg.compute != Config::CPU) {
        std::string reason;
        if (cuda_available(&reason)) {
            try {
                gpu = make_cuda_focus(observed, n, cfg);
            } catch (const std::exception &error) {
                if (cfg.gpu_strict)
                    throw;
                log(cfg, std::string("S3 CUDA focus unavailable; CPU fallback: ") + error.what());
            }
        } else if (cfg.gpu_strict)
            throw std::runtime_error(reason);
    }
    log(cfg, "S3 subset-local duplication scope: " + std::to_string(disco.trees.size()) +
                 " observed components retained; " + (gpu ? "CUDA" : "CPU") + " placements");
    JavaRandom random(cfg.seed);
    std::vector<Tree> proposals;
    proposals.reserve(cfg.rounds);
    for (int round = 0; round < cfg.rounds; ++round) {
        std::vector<Id> order(n);
        std::iota(order.begin(), order.end(), 0);
        random.shuffle(order);
        size_t m = std::max<size_t>(
            4, std::min(n, static_cast<size_t>(std::sqrt(double(n)) * std::log(double(n)) /
                                               std::log(2.0) / 4.0)));
        std::vector<Id> selected(order.begin(), order.begin() + m);
        GreedyTree tree(n);
        {
            Registry subset;
            for (Id taxon : selected)
                subset.add(registry.names[taxon]);
            std::vector<Tree> restricted_genes = restrict_genes(genes, selected, n);
            std::vector<Tree> guides = backbone_guides(disco.trees, selected, n, cfg);
            CandidateSpace candidates(m);
            for (const Tree &guide : guides)
                candidates.add_tree(guide, true);
            if (cfg.cross_tree)
                candidates.add_cross(cfg);
            candidates.finalize();
            Observations backbone_observed = build_observations(restricted_genes);
            InferenceResult inferred = infer(candidates, backbone_observed, subset, cfg);
            tree.import(inferred.tree, selected);
        } // Release backbone trees, cluster tables, transitions and scoring buffers.
        std::vector<Id> edges = tree.preorder;
        std::vector<Id> edge_id(tree.nodes.size(), none);
        std::vector<std::vector<Id>> groups(edges.size());
        for (Id i = 0; i < edges.size(); ++i)
            edge_id[edges[i]] = i;
        const std::vector<Id> backbone_taxa = tree.leaves;
        for (size_t at = m; at < n; ++at) {
            Id taxon = order[at];
            Focus focused = gpu ? gpu->focus(taxon, backbone_taxa)
                                : cpu_focus(observed, n, taxon, backbone_taxa);
            Id edge = locate(tree, taxon, focused, random);
            groups[edge_id[edge]].push_back(taxon);
        }
        std::vector<Matrix> matrices = partition_matrices(disco.trees, groups, n, cfg);
        size_t nonempty = 0;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (groups[i].empty())
                continue;
            ++nonempty;
            Id edge = edges[i], parent = tree.nodes[edge].parent;
            Id partition = group_upgma(tree, groups[i], matrices[i]);
            matrices[i] = Matrix();
            Id joined = tree.joint(edge, partition);
            tree.nodes[joined].parent = parent;
            if (parent == none)
                tree.root = joined;
            else if (tree.nodes[parent].left == edge)
                tree.nodes[parent].left = joined;
            else if (tree.nodes[parent].right == edge)
                tree.nodes[parent].right = joined;
            else
                throw std::runtime_error("S3 backbone parent changed during grouped attachment");
        }
        tree.refresh();
        log(cfg, "S3 round " + std::to_string(round + 1) + "/" + std::to_string(cfg.rounds) + ": " +
                     std::to_string(n - m) + " placements in " + std::to_string(nonempty) +
                     " groups");
        refine_nni(tree, colour, cfg);
        proposals.push_back(tree.export_tree());
    }
    return proposals;
}

} // namespace stelar
