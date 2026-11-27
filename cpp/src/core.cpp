#include "stelar/core.hpp"
#include <cctype>
#include <condition_variable>
#include <exception>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

namespace stelar {
void trim(Bits &a) {
    while (!a.empty() && !a.back())
        a.pop_back();
}
Bits singleton(Id x) {
    Bits a(x / 64 + 1);
    a[x / 64] = 1ULL << (x % 64);
    return a;
}
Bits full_set(size_t n) {
    Bits a((n + 63) / 64, ~0ULL);
    if (n % 64)
        a.back() = (1ULL << (n % 64)) - 1;
    return a;
}
Bits set_union(const Bits &a, const Bits &b) {
    Bits c(std::max(a.size(), b.size()));
    for (size_t i = 0; i < c.size(); ++i)
        c[i] = (i < a.size() ? a[i] : 0) | (i < b.size() ? b[i] : 0);
    return c;
}
Bits set_difference(const Bits &a, const Bits &b) {
    Bits c = a;
    for (size_t i = 0; i < std::min(c.size(), b.size()); ++i)
        c[i] &= ~b[i];
    trim(c);
    return c;
}
Bits set_intersection(const Bits &a, const Bits &b) {
    Bits c(std::min(a.size(), b.size()));
    for (size_t i = 0; i < c.size(); ++i)
        c[i] = a[i] & b[i];
    trim(c);
    return c;
}
bool contains(const Bits &a, Id x) {
    return x / 64 < a.size() && ((a[x / 64] >> (x % 64)) & 1);
}
bool subset(const Bits &a, const Bits &b) {
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i] & ~(i < b.size() ? b[i] : 0))
            return false;
    return true;
}
uint32_t cardinality(const Bits &a) {
    uint32_t c = 0;
    for (auto w : a)
        c += __builtin_popcountll(w);
    return c;
}
uint32_t intersection_size(const Bits &a, const Bits &b) {
    uint32_t c = 0;
    for (size_t i = 0; i < std::min(a.size(), b.size()); ++i)
        c += __builtin_popcountll(a[i] & b[i]);
    return c;
}
std::vector<Id> members(const Bits &a) {
    std::vector<Id> out;
    out.reserve(cardinality(a));
    for (size_t i = 0; i < a.size(); ++i) {
        uint64_t w = a[i];
        while (w) {
            out.push_back(Id(i * 64 + __builtin_ctzll(w)));
            w &= w - 1;
        }
    }
    return out;
}
std::string score_string(Score a) {
    if (!a)
        return "0";
    std::string s;
    while (a) {
        s.push_back('0' + a % 10);
        a /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}
Id Registry::add(const std::string &s) {
    auto it = ids.find(s);
    if (it != ids.end())
        return it->second;
    if (names.size() >= none)
        throw std::length_error("Too many species");
    Id id = names.size();
    names.push_back(s);
    ids.emplace(s, id);
    return id;
}
Id Registry::find(const std::string &s) const {
    auto it = ids.find(s);
    return it == ids.end() ? none : it->second;
}
Id Tree::add_leaf(Id t) {
    if (nodes.size() >= none)
        throw std::length_error("Too many nodes");
    Node a;
    a.taxon = t;
    nodes.push_back(a);
    return nodes.size() - 1;
}
Id Tree::add_internal(const std::vector<Id> &c, Event e) {
    if (c.empty())
        throw std::invalid_argument("Internal node without children");
    if (nodes.size() >= none)
        throw std::length_error("Too many nodes");
    Node a;
    a.child = c.front();
    a.event = e;
    Id id = nodes.size();
    nodes.push_back(a);
    for (size_t i = 0; i < c.size(); ++i) {
        nodes[c[i]].parent = id;
        nodes[c[i]].next = i + 1 < c.size() ? c[i + 1] : none;
    }
    return id;
}
std::vector<Id> Tree::children(Id i) const {
    std::vector<Id> c;
    for (Id j = nodes.at(i).child; j != none; j = nodes[j].next)
        c.push_back(j);
    return c;
}
std::vector<Id> Tree::postorder() const {
    if (root == none)
        return {};
    std::vector<Id> a;
    std::vector<std::pair<Id, bool>> s{{root, false}};
    while (!s.empty()) {
        auto [i, done] = s.back();
        s.pop_back();
        if (done) {
            a.push_back(i);
            continue;
        }
        s.push_back({i, true});
        auto ch = children(i);
        for (auto j = ch.rbegin(); j != ch.rend(); ++j)
            s.push_back({*j, false});
    }
    return a;
}
void Tree::index() {
    leaves.clear();
    if (root == none)
        return;
    std::vector<std::pair<Id, bool>> s{{root, false}};
    while (!s.empty()) {
        auto [i, done] = s.back();
        s.pop_back();
        if (!done) {
            nodes[i].lo = leaves.size();
            if (nodes[i].child == none) {
                leaves.push_back(nodes[i].taxon);
                nodes[i].hi = leaves.size();
            } else {
                s.push_back({i, true});
                auto c = children(i);
                for (auto j = c.rbegin(); j != c.rend(); ++j)
                    s.push_back({*j, false});
            }
        } else
            nodes[i].hi = leaves.size();
    }
}
void Tree::compact() {
    if (root == none) {
        nodes.clear();
        leaves.clear();
        return;
    }
    Tree t;
    std::vector<Id> map(nodes.size(), none);
    for (Id i : postorder()) {
        if (nodes[i].child == none)
            map[i] = t.add_leaf(nodes[i].taxon);
        else {
            std::vector<Id> c;
            for (Id j : children(i))
                c.push_back(map[j]);
            map[i] = t.add_internal(c, nodes[i].event);
        }
    }
    t.root = map[root];
    t.missing = std::move(missing);
    t.scoped = scoped;
    t.index();
    *this = std::move(t);
}
bool Tree::binary() const {
    for (Id i : postorder())
        if (nodes[i].child != none && children(i).size() != 2)
            return false;
    return true;
}
Bits Tree::taxa() const {
    Bits b;
    for (Id x : leaves) {
        if (b.size() <= x / 64)
            b.resize(x / 64 + 1);
        b[x / 64] |= 1ULL << (x % 64);
    }
    return b;
}

namespace {
struct Parser {
    const std::string &s;
    Registry &r;
    bool reg;
    size_t p = 0;
    Tree tree;
    void ws() {
        while (p < s.size()) {
            if (std::isspace(static_cast<unsigned char>(s[p]))) {
                ++p;
                continue;
            }
            if (s[p] == '[') {
                int d = 1;
                ++p;
                while (p < s.size() && d) {
                    if (s[p] == '[')
                        ++d;
                    if (s[p] == ']')
                        --d;
                    ++p;
                }
                if (d)
                    throw std::invalid_argument("Unclosed Newick comment");
                continue;
            }
            break;
        }
    }
    std::string label() {
        ws();
        std::string a;
        if (p < s.size() && s[p] == '\'') {
            ++p;
            while (p < s.size()) {
                char c = s[p++];
                if (c == '\'') {
                    if (p < s.size() && s[p] == '\'') {
                        a += '\'';
                        ++p;
                    } else
                        return a;
                } else
                    a += c;
            }
            throw std::invalid_argument("Unclosed Newick quote");
        }
        while (p < s.size() && s[p] != '(' && s[p] != ')' && s[p] != ',' && s[p] != ':' &&
               s[p] != ';' && s[p] != '[') {
            a += s[p++];
        }
        while (!a.empty() && std::isspace(static_cast<unsigned char>(a.back())))
            a.pop_back();
        return a;
    }
    void length() {
        ws();
        if (p < s.size() && s[p] == ':') {
            ++p;
            while (p < s.size() && s[p] != ',' && s[p] != ')' && s[p] != ';' && s[p] != '[')
                ++p;
            ws();
        }
    }
    Tree parse() {
        std::vector<std::vector<Id>> frames;
        Id completed = none;
        bool expect = true;
        while (true) {
            ws();
            if (p >= s.size())
                break;
            char c = s[p];
            if (expect) {
                if (c == '(') {
                    ++p;
                    frames.emplace_back();
                    continue;
                }
                auto a = label();
                if (a.empty())
                    throw std::invalid_argument("Empty Newick leaf");
                Id tax = reg ? r.add(a) : r.find(a);
                if (tax == none)
                    throw std::invalid_argument("Unknown species: " + a);
                completed = tree.add_leaf(tax);
                length();
                expect = false;
            } else if (c == ',') {
                if (frames.empty())
                    throw std::invalid_argument("Unexpected Newick comma");
                frames.back().push_back(completed);
                completed = none;
                ++p;
                expect = true;
            } else if (c == ')') {
                if (frames.empty())
                    throw std::invalid_argument("Unexpected Newick close");
                frames.back().push_back(completed);
                auto children = std::move(frames.back());
                frames.pop_back();
                ++p;
                auto tag = label();
                completed =
                    children.size() == 1
                        ? children[0]
                        : tree.add_internal(children, tag == "D" ? DUPLICATION : SPECIATION);
                length();
            } else if (c == ';') {
                ++p;
                break;
            } else
                throw std::invalid_argument("Unexpected Newick input at character " +
                                            std::to_string(p));
        }
        if (expect || !frames.empty() || completed == none)
            throw std::invalid_argument("Incomplete Newick tree");
        ws();
        if (p != s.size())
            throw std::invalid_argument("Trailing content after Newick tree");
        tree.root = completed;
        tree.nodes[completed].parent = none;
        tree.index();
        return std::move(tree);
    }
};
uint64_t mix(uint64_t x) {
    x ^= x >> 30;
    x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27;
    x *= 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
uint64_t words_hash(const uint64_t *b, size_t n) {
    uint64_t h = mix(n);
    for (size_t i = 0; i < n; ++i)
        h = mix(h ^ mix(b[i]));
    return h;
}
uint64_t first_hash(const SetPool &s, Id i) {
    uint64_t h = 0, seed = mix(0xDEADBEEFCAFEULL);
    for (Id t : members(s.bits(i))) {
        uint64_t v = mix(t + seed);
        h += v ? v : 1;
    }
    return h;
}
std::string quote(const std::string &s) {
    if (s.find_first_of(" ()[],;:'\t\r\n") == std::string::npos)
        return s;
    std::string a = "'";
    for (char c : s) {
        a += c;
        if (c == '\'')
            a += '\'';
    }
    return a + "'";
}
} // namespace
Tree parse_newick(const std::string &s, Registry &r, bool reg, bool) {
    return Parser{s, r, reg, 0, Tree{}}.parse();
}
std::vector<Tree> read_trees(const std::string &path, Registry &r, bool reg, bool strip) {
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Cannot open tree file: " + path);
    std::vector<Tree> out;
    std::string s;
    char c;
    bool quoted = false;
    int comment = 0;
    while (f.get(c)) {
        s += c;
        if (c == '\'' && !comment) {
            if (quoted && f.peek() == '\'') {
                s += char(f.get());
                continue;
            }
            quoted = !quoted;
        }
        if (!quoted) {
            if (c == '[')
                ++comment;
            else if (c == ']')
                --comment;
            if (c == ';' && !comment) {
                out.push_back(parse_newick(s, r, reg, strip));
                s.clear();
            }
        }
    }
    if (s.find_first_not_of(" \r\n\t") != std::string::npos)
        out.push_back(parse_newick(s, r, reg, strip));
    if (out.empty())
        throw std::invalid_argument("Empty gene tree file");
    return out;
}
std::string newick(const Tree &t, const Registry &r, bool tags) {
    if (t.root == none)
        return ";";
    struct Item {
        Id id;
        char text;
    };
    std::vector<Item> s{{t.root, 0}};
    std::string out;
    while (!s.empty()) {
        auto a = s.back();
        s.pop_back();
        if (a.text) {
            out += a.text;
            continue;
        }
        auto &node = t.nodes[a.id];
        if (node.child == none) {
            out += quote(r.names.at(node.taxon));
            continue;
        }
        if (tags && node.event == DUPLICATION)
            s.push_back({none, 'D'});
        s.push_back({none, ')'});
        auto ch = t.children(a.id);
        for (size_t i = ch.size(); i--;) {
            s.push_back({ch[i], 0});
            if (i)
                s.push_back({none, ','});
        }
        s.push_back({none, '('});
    }
    return out + ";";
}
Tree restrict_tree(const Tree &t, const std::vector<int> &map) {
    Tree o;
    std::vector<Id> m(t.nodes.size(), none);
    for (Id i : t.postorder()) {
        auto &v = t.nodes[i];
        if (v.child == none) {
            if (v.taxon < map.size() && map[v.taxon] >= 0)
                m[i] = o.add_leaf(map[v.taxon]);
        } else {
            std::vector<Id> c;
            for (Id j : t.children(i))
                if (m[j] != none)
                    c.push_back(m[j]);
            if (c.size() == 1)
                m[i] = c[0];
            else if (!c.empty())
                m[i] = o.add_internal(c, v.event);
        }
    }
    if (t.root != none)
        o.root = m[t.root];
    if (o.root != none) {
        o.nodes[o.root].parent = none;
        o.nodes[o.root].next = none;
    }
    o.scoped = t.scoped;
    for (Id x : members(t.missing))
        if (x < map.size() && map[x] >= 0)
            o.missing = set_union(o.missing, singleton(map[x]));
    o.index();
    o.compact();
    return o;
}
std::vector<Tree> restrict_genes(const std::vector<Tree> &g, const std::vector<Id> &sel, size_t n) {
    std::vector<int> map(n, -1);
    for (size_t i = 0; i < sel.size(); ++i)
        map.at(sel[i]) = i;
    std::vector<Tree> o;
    for (auto &t : g) {
        auto a = restrict_tree(t, map);
        if (a.root != none)
            o.push_back(std::move(a));
    }
    return o;
}
void log(const Config &c, const std::string &s) {
    if (!c.quiet && c.verbosity > 0)
        std::cerr << s << '\n';
}
namespace {
thread_local bool in_pool_job = false;
class WorkerPool {
    std::mutex invocation, state;
    std::condition_variable ready, finished;
    std::vector<std::thread> workers;
    std::function<void(size_t)> function;
    std::atomic<size_t> next{0};
    std::atomic<bool> cancelled{false};
    std::exception_ptr error;
    size_t count = 0, pending = 0, generation = 0;
    bool stopping = false;
    void consume() {
        bool previous = in_pool_job;
        in_pool_job = true;
        try {
            while (!cancelled) {
                size_t i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= count)
                    break;
                function(i);
            }
        } catch (...) {
            cancelled = true;
            std::lock_guard<std::mutex> lock(state);
            if (!error)
                error = std::current_exception();
        }
        in_pool_job = previous;
    }
    void stop() {
        {
            std::lock_guard<std::mutex> lock(state);
            stopping = true;
        }
        ready.notify_all();
        for (auto &thread : workers)
            thread.join();
        workers.clear();
    }

  public:
    ~WorkerPool() {
        stop();
    }
    void run(size_t n, int threads, const std::function<void(size_t)> &fn) {
        std::unique_lock<std::mutex> one_invocation(invocation);
        if (workers.size() != size_t(threads - 1)) {
            stop();
            stopping = false;
            generation = 0;
            for (int i = 1; i < threads; ++i)
                workers.emplace_back([this] {
                    size_t seen = 0;
                    while (true) {
                        {
                            std::unique_lock<std::mutex> lock(state);
                            ready.wait(lock, [&] { return stopping || generation != seen; });
                            if (stopping)
                                return;
                            seen = generation;
                        }
                        consume();
                        {
                            std::lock_guard<std::mutex> lock(state);
                            if (!--pending)
                                finished.notify_one();
                        }
                    }
                });
        }
        {
            std::lock_guard<std::mutex> lock(state);
            function = fn;
            count = n;
            next = 0;
            cancelled = false;
            error = nullptr;
            pending = workers.size();
            ++generation;
        }
        ready.notify_all();
        consume();
        {
            std::unique_lock<std::mutex> lock(state);
            finished.wait(lock, [&] { return pending == 0; });
            function = {};
            if (error)
                std::rethrow_exception(error);
        }
    }
};
} // namespace
void parallel_for(size_t n, int threads, const std::function<void(size_t)> &fn) {
    if (!n)
        return;
    threads = std::max(1, std::min(threads, int(std::min<size_t>(n, INT32_MAX))));
    if (threads == 1 || in_pool_job) {
        for (size_t i = 0; i < n; ++i)
            fn(i);
        return;
    }
    static WorkerPool pool;
    pool.run(n, threads, fn);
}
uint32_t JavaRandom::next(int b) {
    state = (state * 0x5DEECE66DULL + 0xB) & ((1ULL << 48) - 1);
    return state >> (48 - b);
}
int JavaRandom::next_int(int bound) {
    if (bound <= 0)
        throw std::invalid_argument("Random bound");
    if ((bound & -bound) == bound)
        return int((uint64_t(bound) * next(31)) >> 31);
    uint32_t b, v;
    do {
        b = next(31);
        v = b % bound;
    } while (b - v + uint32_t(bound - 1) >= 0x80000000U);
    return v;
}
int64_t JavaRandom::next_long() {
    int64_t hi = int32_t(next(32)), lo = int32_t(next(32));
    return int64_t((uint64_t(hi) << 32) + uint64_t(lo));
}
void JavaRandom::shuffle(std::vector<Id> &v) {
    for (size_t i = v.size(); i > 1; --i)
        std::swap(v[i - 1], v[next_int(i)]);
}
void SetPool::rehash(size_t capacity) {
    std::vector<Id> fresh(capacity, none);
    for (Id id = 0; id < records.size(); ++id) {
        auto &r = records[id];
        size_t p = words_hash(arena.data() + r.offset, r.words) & (capacity - 1);
        while (fresh[p] != none)
            p = (p + 1) & (capacity - 1);
        fresh[p] = id;
    }
    slots.swap(fresh);
}
Id SetPool::find(const Bits &raw) const {
    if (slots.empty())
        return none;
    size_t n = raw.size();
    while (n && !raw[n - 1])
        --n;
    size_t p = words_hash(raw.data(), n) & (slots.size() - 1);
    while (slots[p] != none) {
        Id i = slots[p];
        auto &r = records[i];
        if (r.words == n && std::equal(raw.begin(), raw.begin() + n, arena.begin() + r.offset))
            return i;
        p = (p + 1) & (slots.size() - 1);
    }
    return none;
}
Id SetPool::intern(const Bits &raw) {
    if (slots.empty() || (records.size() + 1) * 10 > slots.size() * 7)
        rehash(slots.empty() ? 16 : slots.size() * 2);
    size_t n = raw.size();
    while (n && !raw[n - 1])
        --n;
    size_t p = words_hash(raw.data(), n) & (slots.size() - 1);
    while (slots[p] != none) {
        Id i = slots[p];
        auto &r = records[i];
        if (r.words == n && std::equal(raw.begin(), raw.begin() + n, arena.begin() + r.offset))
            return i;
        p = (p + 1) & (slots.size() - 1);
    }
    if (records.size() >= none)
        throw std::length_error("Too many unique species sets");
    uint32_t count = 0;
    for (size_t j = 0; j < n; ++j)
        count += __builtin_popcountll(raw[j]);
    Id id = records.size();
    records.push_back({arena.size(), uint32_t(n), count});
    arena.insert(arena.end(), raw.begin(), raw.begin() + n);
    slots[p] = id;
    return id;
}
Bits SetPool::bits(Id i) const {
    auto &r = records.at(i);
    return Bits(arena.begin() + r.offset, arena.begin() + r.offset + r.words);
}
uint32_t SetPool::intersect(Id a, const SetPool &b, Id c) const {
    auto &r = records[a];
    auto &s = b.records[c];
    uint32_t n = 0;
    for (size_t i = 0; i < std::min(r.words, s.words); ++i)
        n += __builtin_popcountll(arena[r.offset + i] & b.arena[s.offset + i]);
    return n;
}
CandidateSpace::CandidateSpace(size_t n) : root(none), taxa(n) {
    root = add_cluster(full_set(n), false);
}
Id CandidateSpace::add_cluster(const Bits &b, bool candidate) {
    Id id = sets.intern(b);
    eligible.resize(sets.records.size(), false);
    if (candidate && id != root && sets.count(id))
        eligible[id] = true;
    return id;
}
void CandidateSpace::add_split(Id p, Id a, Id b) {
    if (!sets.count(a) || !sets.count(b) || sets.count(a) + sets.count(b) != sets.count(p))
        return;
    if (a > b)
        std::swap(a, b);
    uint64_t key = (uint64_t(a) << 32) | b;
    if (split_keys.empty() || (split_key_count + 1) * 10 > split_keys.size() * 7) {
        std::vector<uint64_t> fresh(split_keys.empty() ? 16 : split_keys.size() * 2, ~0ULL);
        for (auto k : split_keys)
            if (k != ~0ULL) {
                size_t q = mix(k) & (fresh.size() - 1);
                while (fresh[q] != ~0ULL)
                    q = (q + 1) & (fresh.size() - 1);
                fresh[q] = k;
            }
        split_keys.swap(fresh);
    }
    size_t q = mix(key) & (split_keys.size() - 1);
    while (split_keys[q] != ~0ULL) {
        if (split_keys[q] == key)
            return;
        q = (q + 1) & (split_keys.size() - 1);
    }
    split_keys[q] = key;
    ++split_key_count;
    if (first_hash(sets, b) < first_hash(sets, a))
        std::swap(a, b);
    if (splits.size() >= none)
        throw std::length_error("Too many candidate splits");
    splits.push_back({p, a, b});
    offsets.clear();
    transitions.clear();
}
void CandidateSpace::add_tree(const Tree &t, bool original) {
    std::vector<Bits> live(t.nodes.size());
    for (Id i : t.postorder()) {
        const auto &v = t.nodes[i];
        Bits u;
        if (v.child == none)
            u = singleton(v.taxon);
        else
            for (Id j : t.children(i))
                u = set_union(u, live[j]);
        bool admitted = v.child == none || !original || v.event == SPECIATION;
        Id p = none;
        if (admitted)
            p = add_cluster(u, i != t.root || t.leaves.size() < taxa);
        auto ch = t.children(i);
        if (admitted && ch.size() == 2) {
            Id a = add_cluster(live[ch[0]]), b = add_cluster(live[ch[1]]);
            add_split(p, a, b);
        } else if (admitted && ch.size() > 2) {
            for (Id j : ch)
                if (t.nodes[j].child != none && t.nodes[j].event != SPECIATION)
                    add_cluster(live[j]);
        }
        for (Id j : ch)
            Bits().swap(live[j]);
        live[i] = std::move(u);
    }
}
void CandidateSpace::add_cross(const Config &cfg) {
    if (cfg.compute == Config::GPU) {
        try {
            auto additional = cuda_cross_candidates(*this, cfg);
            for (const auto &split : additional)
                add_split(split.parent, split.left, split.right);
            return;
        } catch (const std::exception &error) {
            if (cfg.gpu_strict)
                throw;
            log(cfg, std::string("CUDA cross-tree fallback: ") + error.what());
        }
    }
    std::vector<std::vector<Id>> bins(taxa + 1);
    std::vector<Id> parents;
    for (Id i = 0; i < eligible.size(); ++i)
        if (eligible[i]) {
            bins[sets.count(i)].push_back(i);
            parents.push_back(i);
        }
    parents.push_back(root);
    std::vector<std::vector<std::pair<Id, Id>>> out(parents.size());
    parallel_for(parents.size(), cfg.threads, [&](size_t k) {
        Id p = parents[k];
        Bits a = sets.bits(p);
        for (size_t s = 1; s <= sets.count(p) / 2; ++s)
            for (Id b : bins[s]) {
                Bits bb = sets.bits(b);
                if (!subset(bb, a))
                    continue;
                Id r = sets.find(set_difference(a, bb));
                if (r != none && eligible[r])
                    out[k].emplace_back(b, r);
            }
    });
    for (size_t k = 0; k < parents.size(); ++k)
        for (auto [a, b] : out[k])
            add_split(parents[k], a, b);
}
void CandidateSpace::finalize() {
    offsets.assign(sets.records.size() + 1, 0);
    for (auto &s : splits)
        ++offsets[s.parent + 1];
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    transitions.resize(splits.size());
    auto cursor = offsets;
    for (Id i = 0; i < splits.size(); ++i)
        transitions[cursor[splits[i].parent]++] = i;
}
std::vector<Id> CandidateSpace::reachable(bool prune) const {
    std::vector<Id> out;
    if (!prune) {
        out.resize(splits.size());
        std::iota(out.begin(), out.end(), 0);
        return out;
    }
    std::vector<uint8_t> seen(sets.records.size());
    std::vector<Id> todo{root};
    seen[root] = 1;
    while (!todo.empty()) {
        Id p = todo.back();
        todo.pop_back();
        for (size_t k = offsets[p]; k < offsets[p + 1]; ++k) {
            Id si = transitions[k];
            out.push_back(si);
            for (Id c : {splits[si].left, splits[si].right})
                if (!seen[c]) {
                    seen[c] = 1;
                    todo.push_back(c);
                }
        }
    }
    return out;
}
Observations build_observations(const std::vector<Tree> &trees) {
    Observations o;
    std::unordered_map<uint64_t, std::vector<Id>> index;
    for (auto &t : trees) {
        std::vector<Bits> live(t.nodes.size());
        Bits all = t.taxa();
        for (Id i : t.postorder()) {
            const auto &v = t.nodes[i];
            if (v.child == none) {
                live[i] = singleton(v.taxon);
                continue;
            }
            auto ch = t.children(i);
            Bits u;
            for (Id j : ch)
                u = set_union(u, live[j]);
            if (v.event == SPECIATION && ch.size() >= 2 &&
                !(ch.size() > 2 && set_difference(all, u).empty())) {
                std::vector<Id> ids;
                for (Id j : ch)
                    ids.push_back(o.sets.intern(live[j]));
                std::sort(ids.begin(), ids.end());
                uint64_t h = mix(ids.size());
                for (Id a : ids)
                    h = mix(h ^ a);
                auto &bucket = index[h];
                Id known = none;
                for (Id p : bucket) {
                    auto &r = o.partitions[p];
                    if (r.children == ids.size() &&
                        std::equal(ids.begin(), ids.end(), o.child_sets.begin() + r.offset)) {
                        known = p;
                        break;
                    }
                }
                if (known != none) {
                    if (o.partitions[known].frequency == UINT32_MAX)
                        throw std::overflow_error("Partition frequency exceeds uint32");
                    ++o.partitions[known].frequency;
                } else {
                    bucket.push_back(o.partitions.size());
                    o.partitions.push_back({o.child_sets.size(), uint32_t(ids.size()), 1});
                    o.child_sets.insert(o.child_sets.end(), ids.begin(), ids.end());
                }
            }
            for (Id j : ch)
                Bits().swap(live[j]);
            live[i] = std::move(u);
        }
    }
    return o;
}
Score split_weight(const SetPool &s, Id a, Id b, const Observations &o) {
    Score two = 0;
    for (auto &p : o.partitions) {
        Score val = 0;
        if (p.children == 2) {
            Id c = o.child_sets[p.offset], d = o.child_sets[p.offset + 1];
            uint64_t ac = s.intersect(a, o.sets, c), ad = s.intersect(a, o.sets, d),
                     bc = s.intersect(b, o.sets, c), bd = s.intersect(b, o.sets, d);
            if (ac && bd)
                val += Score(ac) * bd * (ac + bd - 2);
            if (ad && bc)
                val += Score(ad) * bc * (ad + bc - 2);
        } else {
            uint64_t sa = 0, sb = 0;
            for (size_t j = 0; j < p.children; ++j) {
                Id c = o.child_sets[p.offset + j];
                sa += s.intersect(a, o.sets, c);
                sb += s.intersect(b, o.sets, c);
            }
            for (size_t j = 0; j < p.children; ++j) {
                Id c = o.child_sets[p.offset + j];
                uint64_t ac = s.intersect(a, o.sets, c), bc = s.intersect(b, o.sets, c);
                if (ac >= 2)
                    val += Score(ac) * (ac - 1) * (sb - bc);
                if (bc >= 2)
                    val += Score(bc) * (bc - 1) * (sa - ac);
            }
        }
        two += val * p.frequency;
    }
    return two / 2;
}
InferenceResult infer(CandidateSpace &c, const Observations &o, const Registry &r,
                      const Config &cfg) {
    c.finalize();
    auto ids = c.reachable(cfg.prune);
    auto start = std::chrono::steady_clock::now();
    std::vector<Score> w(c.splits.size());
    bool gpu = cfg.compute == Config::GPU;
    std::vector<Score> weights;
    if (gpu) {
        try {
            weights = cuda_weights(c, o, ids, cfg);
        } catch (const std::exception &e) {
            if (cfg.gpu_strict)
                throw;
            log(cfg, std::string("CUDA weights failed; CPU fallback: ") + e.what());
            gpu = false;
        }
    }
    if (!gpu) {
        weights.resize(ids.size());
        parallel_for(ids.size(), cfg.threads, [&](size_t k) {
            const auto &s = c.splits[ids[k]];
            weights[k] = split_weight(c.sets, s.left, s.right, o);
        });
    }
    for (size_t k = 0; k < ids.size(); ++k)
        w[ids[k]] = weights[k];
    log(cfg,
        "Weight calculation: " + std::to_string(ids.size()) + " scored splits; " +
            (gpu ? "CUDA" : "CPU") + " " +
            std::to_string(
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()) +
            " s");
    if (cfg.option("--verify-weights")) {
        std::string path = cfg.output + ".weights.tsv";
        std::ofstream f(path);
        if (!f)
            throw std::runtime_error("Cannot open weight diagnostic: " + path);
        for (Id i : ids) {
            auto &s = c.splits[i];
            f << i << '\t' << score_string(w[i]) << '\t';
            for (Id t : members(c.sets.bits(s.left)))
                f << r.names[t] << ',';
            f << '\t';
            for (Id t : members(c.sets.bits(s.right)))
                f << r.names[t] << ',';
            f << '\n';
        }
    }
    std::vector<Id> order(c.sets.records.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](Id a, Id b) { return c.sets.count(a) < c.sets.count(b); });
    std::vector<Score> dp(order.size());
    std::vector<Id> chosen(order.size(), none);
    for (Id p : order) {
        bool found = false;
        for (size_t k = c.offsets[p]; k < c.offsets[p + 1]; ++k) {
            Id i = c.transitions[k];
            auto &s = c.splits[i];
            Score v = w[i] + dp[s.left] + dp[s.right];
            if (!found || v > dp[p]) {
                found = true;
                dp[p] = v;
                chosen[p] = i;
            }
        }
    }
    Tree t;
    std::vector<std::pair<Id, bool>> stack{{c.root, false}};
    std::vector<Id> node(order.size(), none);
    while (!stack.empty()) {
        auto [p, done] = stack.back();
        stack.pop_back();
        Id i = chosen[p];
        if (i == none) {
            std::vector<Id> ch;
            for (Id tax : members(c.sets.bits(p)))
                ch.push_back(t.add_leaf(tax));
            node[p] = ch.size() == 1 ? ch[0] : t.add_internal(ch);
            continue;
        }
        auto &s = c.splits[i];
        if (!done) {
            stack.push_back({p, true});
            stack.push_back({s.right, false});
            stack.push_back({s.left, false});
        } else
            node[p] = t.add_internal({node[s.left], node[s.right]});
    }
    t.root = node[c.root];
    t.index();
    size_t n = 0;
    for (auto e : c.eligible)
        n += e;
    return {std::move(t), dp[c.root], n, c.splits.size(), ids.size()};
}
Focus cpu_focus(const Observations &o, size_t n, Id x, const std::vector<Id> &sel) {
    Focus f;
    f.n = n;
    f.selected = sel;
    size_t m = sel.size();
    f.pair.assign(m * m, 0);
    f.sister.assign(m * m, 0);
    std::vector<int> map(n, -1);
    for (size_t i = 0; i < m; ++i)
        map.at(sel[i]) = i;
    for (auto &p : o.partitions)
        for (size_t side = 0; side < p.children; ++side) {
            Id a = o.child_sets[p.offset + side];
            if (!contains(o.sets.bits(a), x))
                continue;
            auto inside = members(o.sets.bits(a));
            for (size_t out = 0; out < p.children; ++out) {
                if (side == out)
                    continue;
                auto outside = members(o.sets.bits(o.child_sets[p.offset + out]));
                for (size_t i = 0; i < outside.size(); ++i) {
                    Id u = outside[i];
                    if (u == x || map[u] < 0)
                        continue;
                    size_t ui = map[u];
                    for (size_t j = i + 1; j < outside.size(); ++j) {
                        Id v = outside[j];
                        if (v != x && map[v] >= 0)
                            f.pair[ui * m + map[v]] += p.frequency;
                    }
                    for (Id v : inside)
                        if (v != x && u != v && map[v] >= 0)
                            f.sister[ui * m + map[v]] += p.frequency;
                }
            }
        }
    return f;
}
uint64_t Focus::support(Id a, Id b, Id c, Id x) const {
    if (a == b || a == c || b == c)
        return 0;
    if (a != x && b != x && c != x)
        throw std::invalid_argument("Triplet excludes focus taxon");
    auto pos = [&](Id t) {
        auto it = std::find(selected.begin(), selected.end(), t);
        if (it == selected.end())
            throw std::out_of_range("Taxon outside focused slice");
        return size_t(it - selected.begin());
    };
    size_t m = selected.size();
    if (c == x) {
        if (a > b)
            std::swap(a, b);
        return pair[pos(a) * m + pos(b)];
    }
    Id other = a == x ? b : a;
    return sister[pos(c) * m + pos(other)];
}
} // namespace stelar
