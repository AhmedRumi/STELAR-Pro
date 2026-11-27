#include "stelar/core.hpp"

#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>

namespace stelar {
namespace {

size_t packed(size_t n, size_t a, size_t b) {
    if (a > b)
        std::swap(a, b);
    return a * n - a * (a + 1) / 2 + b - a - 1;
}

// Temporary integer sets use flat open-addressed tables instead of one heap
// allocation per species. Small bags need only their contiguous member list.
struct TaxaBag {
    std::vector<Id> ids, slots;
    static size_t hash(Id x) {
        uint32_t h = x * 0x9e3779b9U;
        h ^= h >> 16;
        return h;
    }
    void rehash(size_t size) {
        std::vector<Id> replacement(size, none);
        for (Id x : ids) {
            size_t p = hash(x) & (size - 1);
            while (replacement[p] != none)
                p = (p + 1) & (size - 1);
            replacement[p] = x;
        }
        slots.swap(replacement);
    }
    void insert(Id x) {
        if (slots.empty()) {
            if (std::find(ids.begin(), ids.end(), x) != ids.end())
                return;
            ids.push_back(x);
            if (ids.size() == 9)
                rehash(32);
            return;
        }
        if ((ids.size() + 1) * 3 >= slots.size() * 2)
            rehash(slots.size() * 2);
        size_t p = hash(x) & (slots.size() - 1);
        while (slots[p] != none && slots[p] != x)
            p = (p + 1) & (slots.size() - 1);
        if (slots[p] == none) {
            slots[p] = x;
            ids.push_back(x);
        }
    }
};

// Sparse species sets exist only on the active postorder frontier.
// The persistent scope is one trimmed missing-species mask per component.
DiscoResult decompose_family(const Tree &gene, size_t n, Scope scope) {
    if (!gene.binary())
        throw std::invalid_argument("DISCO requires binary rooted/tagged gene trees");
    const auto order = gene.postorder();
    std::vector<TaxaBag> frontier(gene.nodes.size());
    std::vector<Id> detached(gene.nodes.size(), none), roots;
    std::vector<Bits> scopes;
    roots.reserve(order.size() / 2);
    scopes.reserve(order.size() / 2);
    for (Id u : order) {
        const auto &node = gene.nodes[u];
        if (node.child == none) {
            frontier[u].insert(node.taxon);
            continue;
        }
        const auto children = gene.children(u);
        const Id a = children[0], b = children[1];
        auto &left = frontier[a];
        auto &right = frontier[b];
        if (node.event == DUPLICATION) {
            // DISCO cuts the right child when distinct-species counts tie.
            const bool cut_left = left.ids.size() < right.ids.size();
            detached[u] = cut_left ? a : b;
            roots.push_back(detached[u]);
            if (scope != Scope::GLOBAL) {
                Bits saved;
                const auto &source = cut_left ? left : right;
                for (Id x : source.ids) {
                    if (saved.size() <= x / 64)
                        saved.resize(x / 64 + 1);
                    saved[x / 64] |= uint64_t(1) << (x % 64);
                }
                if (scope == Scope::DUPLICATION) {
                    for (Id x : (cut_left ? right : left).ids) {
                        if (saved.size() <= x / 64)
                            saved.resize(x / 64 + 1);
                        saved[x / 64] |= uint64_t(1) << (x % 64);
                    }
                }
                scopes.push_back(std::move(saved));
            }
        }
        if (left.ids.size() < right.ids.size())
            std::swap(left, right);
        for (Id x : right.ids)
            left.insert(x);
        frontier[u] = std::move(left);
        right = TaxaBag{};
    }
    roots.push_back(gene.root);
    if (scope != Scope::GLOBAL) {
        Bits top;
        for (Id x : frontier[gene.root].ids) {
            if (top.size() <= x / 64)
                top.resize(x / 64 + 1);
            top[x / 64] |= uint64_t(1) << (x % 64);
        }
        scopes.push_back(std::move(top));
    }
    std::vector<TaxaBag>().swap(frontier);
    DiscoResult out;
    out.cuts = roots.size() - 1;
    std::vector<Id> translated(gene.nodes.size(), none);
    for (size_t ri = 0; ri < roots.size(); ++ri) {
        Tree component;
        // Iterative reconstruction also suppresses unary cut vertices; the
        // original family and its event tags remain untouched for scoring.
        std::vector<std::pair<Id, bool>> todo{{roots[ri], false}};
        while (!todo.empty()) {
            const auto entry = todo.back();
            todo.pop_back();
            const Id u = entry.first;
            if (gene.nodes[u].child == none) {
                translated[u] = component.add_leaf(gene.nodes[u].taxon);
                continue;
            }
            if (!entry.second) {
                todo.emplace_back(u, true);
                const auto ch = gene.children(u);
                for (auto it = ch.rbegin(); it != ch.rend(); ++it)
                    if (*it != detached[u])
                        todo.emplace_back(*it, false);
            } else {
                std::vector<Id> kept;
                for (Id c : gene.children(u))
                    if (c != detached[u])
                        kept.push_back(translated[c]);
                translated[u] =
                    kept.size() == 1 ? kept[0] : component.add_internal(kept, gene.nodes[u].event);
            }
        }
        component.root = translated[roots[ri]];
        component.nodes[component.root].parent = none;
        component.index();
        if (component.leaves.size() < 4) {
            ++out.discarded;
            continue;
        }
        auto present = component.taxa();
        if (cardinality(present) != component.leaves.size())
            throw std::runtime_error("DISCO component still contains duplicate species");
        component.scoped = scope != Scope::GLOBAL;
        if (component.scoped)
            component.missing = set_difference(scopes[ri], present);
        out.trees.push_back(std::move(component));
    }
    (void)n;
    return out;
}

struct Occurrence {
    uint32_t tree;
    Id node;
};

} // namespace

double Matrix::get(size_t a, size_t b) const {
    if (a >= n || b >= n)
        throw std::out_of_range("similarity species index");
    return a == b ? 1.0 : values[packed(n, a, b)];
}
void Matrix::set(size_t a, size_t b, double value) {
    if (a >= n || b >= n)
        throw std::out_of_range("similarity species index");
    if (a != b)
        values[packed(n, a, b)] = value;
}

DiscoResult decompose(const std::vector<Tree> &genes, size_t n, Scope scope, const Config &cfg) {
    DiscoResult result;
    // Bound live copies/frontier sets while retaining deterministic family and
    // within-family component order regardless of worker scheduling.
    const size_t batch_size = std::max<size_t>(1, std::min<size_t>(32, std::max(1, cfg.threads)));
    for (size_t from = 0; from < genes.size(); from += batch_size) {
        const size_t count = std::min(batch_size, genes.size() - from);
        std::vector<DiscoResult> families(count);
        parallel_for(count, cfg.threads,
                     [&](size_t i) { families[i] = decompose_family(genes[from + i], n, scope); });
        for (auto &family : families) {
            result.cuts += family.cuts;
            result.discarded += family.discarded;
            for (auto &t : family.trees)
                result.trees.push_back(std::move(t));
        }
    }
    log(cfg, "DISCO: " + std::to_string(genes.size()) + " families -> " +
                 std::to_string(result.trees.size()) + " retained components; " +
                 std::to_string(result.cuts) + " duplication cuts; " +
                 std::to_string(result.discarded) + " components below four leaves discarded");
    return result;
}

Matrix triplet_matrix(const std::vector<Tree> &trees, size_t n, const Config &cfg, bool allow_gpu) {
    Matrix result(n);
    if (allow_gpu && cfg.compute != Config::CPU) {
        std::string reason;
        if (cuda_available(&reason)) {
            try {
                if (cuda_triplet_matrix(trees, n, cfg, result))
                    return result;
            } catch (const std::exception &error) {
                if (cfg.gpu_strict)
                    throw;
                log(cfg, std::string("CUDA triplet similarity failed; using CPU: ") + error.what());
                // A failed device call may have partially filled the output.
                std::fill(result.values.begin(), result.values.end(), 0.0);
            }
            if (cfg.gpu_strict)
                throw std::runtime_error("CUDA triplet similarity failed in strict GPU mode");
            log(cfg, "CUDA triplet similarity unavailable for this input; using CPU");
        } else if (cfg.gpu_strict)
            throw std::runtime_error("CUDA triplet similarity unavailable: " + reason);
    }
    std::vector<uint64_t> denominator(result.values.size(), 0);
    // Each species owns its upper-triangle row. Occurrences are sparse CSR;
    // no component allocates an n-species position map or an Euler/RMQ table.
    // For a leaf a, siblings along the a-to-root path partition all other
    // leaves b. Their common parent is exactly LCA(a,b), giving k-|LCA|.
    constexpr size_t occurrence_cap = 1 << 20;
    for (size_t from = 0; from < trees.size();) {
        size_t to = from, total = 0;
        while (to < trees.size() &&
               (to == from || total + trees[to].leaves.size() <= occurrence_cap)) {
            total += trees[to].leaves.size();
            ++to;
        }
        std::vector<size_t> offsets(n + 1, 0);
        for (size_t ti = from; ti < to; ++ti)
            for (Id u : trees[ti].postorder())
                if (trees[ti].nodes[u].child == none)
                    ++offsets[trees[ti].nodes[u].taxon + 1];
        std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
        auto cursor = offsets;
        std::vector<Occurrence> occurrences(total);
        for (size_t ti = from; ti < to; ++ti)
            for (Id u : trees[ti].postorder()) {
                const auto &node = trees[ti].nodes[u];
                if (node.child == none)
                    occurrences[cursor[node.taxon]++] = {uint32_t(ti - from), u};
            }
        parallel_for(n, cfg.threads, [&](size_t a) {
            for (size_t oi = offsets[a]; oi < offsets[a + 1]; ++oi) {
                const auto &occurrence = occurrences[oi];
                const auto &tree = trees[from + occurrence.tree];
                const size_t k = tree.leaves.size();
                if (k < 3)
                    continue;
                Id child = occurrence.node;
                for (Id parent = tree.nodes[child].parent; parent != none;
                     child = parent, parent = tree.nodes[parent].parent) {
                    const auto &p = tree.nodes[parent];
                    const double contribution = double(k - (p.hi - p.lo));
                    for (Id sibling = p.child; sibling != none; sibling = tree.nodes[sibling].next)
                        if (sibling != child) {
                            for (uint32_t pos = tree.nodes[sibling].lo;
                                 pos < tree.nodes[sibling].hi; ++pos) {
                                const Id b = tree.leaves[pos];
                                if (b <= a)
                                    continue;
                                const size_t index = packed(n, a, b);
                                result.values[index] += contribution;
                                denominator[index] += k - 2;
                            }
                        }
                }
            }
        });
        from = to;
    }
    parallel_for(result.values.size(), cfg.threads, [&](size_t i) {
        result.values[i] = denominator[i] ? result.values[i] / double(denominator[i]) : 0.0;
    });
    return result;
}

std::vector<Id> sorted_rows(const Matrix &matrix, const Config &cfg) {
    if (matrix.n && matrix.n > std::numeric_limits<size_t>::max() / matrix.n)
        throw std::length_error("completion anchor table dimensions overflow");
    std::vector<Id> rows(matrix.n * matrix.n);
    parallel_for(matrix.n, cfg.threads, [&](size_t x) {
        auto first = rows.begin() + x * matrix.n;
        std::iota(first, first + matrix.n, Id(0));
        // Java explicitly sorts distance, rather than similarity: subtraction
        // can round distinct near-zero similarities to the same distance.
        std::sort(first, first + matrix.n, [&](Id a, Id b) {
            const double da = 1.0 - matrix.get(x, a), db = 1.0 - matrix.get(x, b);
            return da < db || (da == db && a > b);
        });
    });
    return rows;
}

void complete_tree(Tree &tree, const std::vector<Id> &rows, size_t n, bool scoped) {
    const Bits missing = scoped ? tree.missing : set_difference(full_set(n), tree.taxa());
    if (cardinality(missing) == 0)
        return;
    if (rows.size() != n * n)
        throw std::invalid_argument("completion requires nearest-neighbor rows");
    std::vector<Id> leaf_node(n, none);
    for (Id u : tree.postorder())
        if (tree.nodes[u].child == none)
            leaf_node[tree.nodes[u].taxon] = u;
    tree.nodes.reserve(tree.nodes.size() + 2 * cardinality(missing));
    for (Id x : members(missing)) {
        if (x >= n || leaf_node[x] != none)
            throw std::runtime_error("invalid missing species in DISCO scope");
        Id anchor = none;
        for (size_t rank = 0; rank < n; ++rank) {
            const Id a = rows[size_t(x) * n + rank];
            if (leaf_node[a] != none) {
                anchor = leaf_node[a];
                break;
            }
        }
        if (anchor == none)
            throw std::runtime_error("empty completion component has no anchor");
        const Id old_parent = tree.nodes[anchor].parent, old_next = tree.nodes[anchor].next;
        const Id fresh = tree.add_leaf(x);
        // Java's completion join has no biological event tag. Streaming S4
        // admits it explicitly; S6 backbone hashing still respects this tag.
        const Id join = tree.add_internal({anchor, fresh}, LEAF);
        tree.nodes[join].parent = old_parent;
        tree.nodes[join].next = old_next;
        if (old_parent == none)
            tree.root = join;
        else {
            Id &first = tree.nodes[old_parent].child;
            if (first == anchor)
                first = join;
            else {
                Id prev = first;
                while (prev != none && tree.nodes[prev].next != anchor)
                    prev = tree.nodes[prev].next;
                if (prev == none)
                    throw std::runtime_error("broken DISCO parent link during completion");
                tree.nodes[prev].next = join;
            }
        }
        leaf_node[x] = fresh;
    }
    tree.missing.clear();
    tree.index();
}

Tree upgma(Matrix matrix, const Config &cfg) {
    const size_t n = matrix.n;
    Tree tree;
    if (!n)
        return tree;
    tree.nodes.reserve(2 * n - 1);
    std::vector<Id> active(n), position(n), roots(n), best_neighbor(n, none);
    std::vector<double> best(n, -std::numeric_limits<double>::infinity());
    std::vector<uint32_t> weight(n, 1);
    std::iota(active.begin(), active.end(), Id(0));
    position = active;
    for (Id a = 0; a < n; ++a)
        roots[a] = tree.add_leaf(a);
    auto recompute = [&](Id a) {
        double value = -std::numeric_limits<double>::infinity();
        Id neighbor = none;
        for (Id b : active)
            if (b != a && matrix.get(a, b) > value) {
                value = matrix.get(a, b);
                neighbor = b;
            }
        best[a] = value;
        best_neighbor[a] = neighbor;
    };
    parallel_for(n, n >= 256 ? cfg.threads : 1, [&](size_t a) { recompute(Id(a)); });
    while (active.size() > 1) {
        Id a = none;
        double value = -std::numeric_limits<double>::infinity();
        // Strict comparison, active-array order, and cached-neighbor repairs
        // are deliberate: recomputing all pairs changes Java's equal ties.
        for (Id candidate : active)
            if (best[candidate] > value) {
                value = best[candidate];
                a = candidate;
            }
        if (a == none || best_neighbor[a] == none)
            throw std::runtime_error("UPGMA has no finite pair");
        const Id b = best_neighbor[a];
        roots[a] = tree.add_internal({roots[a], roots[b]});
        const double wa = weight[a], wb = weight[b];
        weight[a] += weight[b];
        const size_t erase_at = position[b];
        const Id last = active.back();
        active[erase_at] = last;
        position[last] = Id(erase_at);
        active.pop_back();
        std::vector<uint8_t> stale(active.size(), 0);
        // These updates are only a handful of operations per active row;
        // avoid launching a worker team for a small UPGMA merge.
        parallel_for(active.size(), active.size() >= 8192 ? cfg.threads : 1, [&](size_t i) {
            const Id k = active[i];
            if (k == a)
                return;
            const double updated = (matrix.get(a, k) * wa + matrix.get(b, k) * wb) / (wa + wb);
            matrix.set(a, k, updated);
            if (best_neighbor[k] == b)
                stale[i] = 1;
            else if (updated > best[k]) {
                best_neighbor[k] = a;
                best[k] = updated;
            } else if (best_neighbor[k] == a && updated < best[k])
                stale[i] = 1;
        });
        const size_t stale_count = std::count(stale.begin(), stale.end(), uint8_t(1));
        parallel_for(active.size(), stale_count * active.size() >= 131072 ? cfg.threads : 1,
                     [&](size_t i) {
                         if (stale[i])
                             recompute(active[i]);
                     });
        recompute(a);
    }
    tree.root = roots[active[0]];
    tree.nodes[tree.root].parent = none;
    tree.index();
    return tree;
}

void add_s4_candidates(CandidateSpace &candidates, const std::vector<Tree> &genes,
                       const Registry &registry, const Config &cfg, Scope requested_scope) {
    if (cfg.exhaustive)
        for (const auto &tree : genes)
            candidates.add_tree(tree, true);
    const Scope scope = cfg.exhaustive ? Scope::GLOBAL : requested_scope;
    auto disco = decompose(genes, registry.size(), scope, cfg);
    if (disco.trees.empty()) {
        if (!cfg.dump_completed.empty()) {
            std::ofstream empty(cfg.dump_completed);
            if (!empty)
                throw std::runtime_error("cannot open completed-tree dump: " + cfg.dump_completed);
        }
        log(cfg, "No DISCO component has four leaves; completion/UPGMA skipped");
        return;
    }
    auto matrix = triplet_matrix(disco.trees, registry.size(), cfg);
    bool needs_completion = false;
    for (const auto &tree : disco.trees)
        if (scope == Scope::GLOBAL ? tree.leaves.size() < registry.size()
                                   : cardinality(tree.missing) > 0) {
            needs_completion = true;
            break;
        }
    auto rows = needs_completion ? sorted_rows(matrix, cfg) : std::vector<Id>{};
    auto guide = upgma(std::move(matrix), cfg);
    std::ofstream dump;
    if (!cfg.dump_completed.empty()) {
        dump.open(cfg.dump_completed);
        if (!dump)
            throw std::runtime_error("cannot open completed-tree dump: " + cfg.dump_completed);
    }
    const size_t batch_size =
        std::max<size_t>(64, std::min<size_t>(256, size_t(std::max(1, cfg.threads)) * 8));
    size_t inserted = 0, excluded = 0, complete = 0;
    const size_t before = candidates.sets.records.size();
    for (size_t from = 0; from < disco.trees.size(); from += batch_size) {
        const size_t count = std::min(batch_size, disco.trees.size() - from);
        for (size_t i = 0; i < count; ++i) {
            const auto &tree = disco.trees[from + i];
            const size_t missing = scope == Scope::GLOBAL ? registry.size() - tree.leaves.size()
                                                          : cardinality(tree.missing);
            inserted += missing;
            excluded += registry.size() - tree.leaves.size() - missing;
            complete += missing == 0;
        }
        parallel_for(count, cfg.threads, [&](size_t i) {
            complete_tree(disco.trees[from + i], rows, registry.size(), scope != Scope::GLOBAL);
        });
        for (size_t i = 0; i < count; ++i) {
            auto &tree = disco.trees[from + i];
            if (dump)
                dump << newick(tree, registry) << '\n';
            candidates.add_tree(tree, false);
            // Sets/transitions own packed species masks. No completed topology
            // or dense exemplar index survives candidate insertion.
            tree = Tree{};
        }
    }
    const size_t after_disco = candidates.sets.records.size();
    candidates.add_tree(guide, false);
    std::ostringstream statistics;
    statistics << "DISCO scopes: " << complete << " components already complete; " << inserted
               << " leaves inserted; " << excluded
               << " outside-scope insertions avoided. Candidates: +" << after_disco - before
               << " DISCO, +" << candidates.sets.records.size() - after_disco << " UPGMA; "
               << candidates.splits.size() << " unique local splits";
    log(cfg, statistics.str());
}

} // namespace stelar
