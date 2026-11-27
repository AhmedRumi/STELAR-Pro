#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
namespace stelar {
using Id = uint32_t;
constexpr Id none = std::numeric_limits<Id>::max();
using Bits = std::vector<uint64_t>;
using Score = unsigned __int128;
void trim(Bits &);
Bits singleton(Id);
Bits full_set(size_t);
Bits set_union(const Bits &, const Bits &);
Bits set_difference(const Bits &, const Bits &);
Bits set_intersection(const Bits &, const Bits &);
bool contains(const Bits &, Id);
bool subset(const Bits &, const Bits &);
uint32_t cardinality(const Bits &);
uint32_t intersection_size(const Bits &, const Bits &);
std::vector<Id> members(const Bits &);
std::string score_string(Score);
struct Registry {
    std::vector<std::string> names;
    std::unordered_map<std::string, Id> ids;
    Id add(const std::string &);
    Id find(const std::string &) const;
    size_t size() const {
        return names.size();
    }
};
enum Event : uint8_t { LEAF = 0, SPECIATION = 1, DUPLICATION = 2 };
struct Node {
    Id parent = none, child = none, next = none, taxon = none;
    uint32_t lo = 0, hi = 0;
    Event event = LEAF;
};
struct Tree {
    std::vector<Node> nodes;
    Id root = none;
    std::vector<Id> leaves;
    Bits missing;
    bool scoped = false;
    Id add_leaf(Id);
    Id add_internal(const std::vector<Id> &, Event = SPECIATION);
    std::vector<Id> children(Id) const;
    std::vector<Id> postorder() const;
    void index();
    void compact();
    bool binary() const;
    Bits taxa() const;
};
Tree parse_newick(const std::string &, Registry &, bool register_names = true,
                  bool strip_copy_suffix = true);
std::vector<Tree> read_trees(const std::string &, Registry &, bool register_names = true,
                             bool strip_copy_suffix = true);
std::string newick(const Tree &, const Registry &, bool event_tags = false);
Tree restrict_tree(const Tree &, const std::vector<int> &);
std::vector<Tree> restrict_genes(const std::vector<Tree> &, const std::vector<Id> &, size_t);
struct Config {
    int mode = 1, threads = 1, seeds = 2, rounds = 4;
    int64_t seed = 0xDEADBEEFCAFE;
    enum Compute { AUTO, CPU, GPU } compute = AUTO;
    bool gpu_strict = false, exhaustive = false, cross_tree = false, prune = true,
         score_only = false, quiet = false;
    int verbosity = 1;
    std::string score_type = "auto";
    std::string input, output, log_file, taxa_file, gene_map, astral_executable, score_tree;
    std::string dump_clusters, dump_completed, completion_method = "similarity";
    bool keep_polytomy = false, autocomplete = false, no_gpu_batch = false;
    int gpu_batch_size = 0, gpu_batches = 0, gpu_tile = 512, gpu_sim_cap_mb = 512,
        gpu_treewalk_cap_mb = 512;
    double gpu_occupancy = .75, gpu_control = 0, gpu_progress = 0;
    int gpu_dp_cap = 10000000, gpu_dp_steps = 8;
    double gpu_dp_progress = 0;
    std::unordered_map<std::string, std::string> options;
    bool option(const std::string &s) const {
        return options.count(s);
    }
};
void log(const Config &, const std::string &);
void parallel_for(size_t, int, const std::function<void(size_t)> &);
struct JavaRandom {
    uint64_t state;
    explicit JavaRandom(int64_t s) : state((uint64_t(s) ^ 0x5DEECE66DULL) & ((1ULL << 48) - 1)) {}
    uint32_t next(int);
    int next_int(int);
    int64_t next_long();
    void shuffle(std::vector<Id> &);
};
class SetPool {
  public:
    struct Record {
        uint64_t offset;
        uint32_t words, count;
    };
    Id intern(const Bits &);
    Id find(const Bits &) const;
    Bits bits(Id) const;
    uint32_t count(Id i) const {
        return records.at(i).count;
    }
    uint32_t intersect(Id, const SetPool &, Id) const;
    const uint64_t *data(Id i) const {
        return arena.data() + records[i].offset;
    }
    std::vector<Record> records;
    std::vector<uint64_t> arena;

  private:
    std::vector<Id> slots;
    void rehash(size_t);
};
struct Split {
    Id parent, left, right;
};
struct CandidateSpace {
    explicit CandidateSpace(size_t);
    SetPool sets;
    Id root;
    size_t taxa;
    std::vector<uint8_t> eligible;
    std::vector<Split> splits;
    std::vector<uint64_t> offsets;
    std::vector<Id> transitions;
    Id add_cluster(const Bits &, bool candidate = true);
    void add_split(Id, Id, Id);
    void add_tree(const Tree &, bool original = false);
    void add_cross(const Config &);
    void finalize();
    std::vector<Id> reachable(bool) const;

  private:
    std::vector<uint64_t> split_keys;
    size_t split_key_count = 0;
};
struct Observation {
    uint64_t offset;
    uint32_t children, frequency;
};
struct Observations {
    SetPool sets;
    std::vector<Id> child_sets;
    std::vector<Observation> partitions;
};
Observations build_observations(const std::vector<Tree> &);
Score split_weight(const SetPool &, Id, Id, const Observations &);
struct InferenceResult {
    Tree tree;
    Score score = 0;
    size_t clusters = 0, splits = 0, scored = 0;
};
InferenceResult infer(CandidateSpace &, const Observations &, const Registry &, const Config &);
enum class Scope { GLOBAL, LINEAGE, DUPLICATION };
struct DiscoResult {
    std::vector<Tree> trees;
    size_t cuts = 0, discarded = 0;
};
DiscoResult decompose(const std::vector<Tree> &, size_t, Scope, const Config &);
struct Matrix {
    size_t n = 0;
    std::vector<double> values;
    explicit Matrix(size_t s = 0) : n(s), values(s * (s ? s - 1 : 0) / 2, 0.0) {}
    double get(size_t, size_t) const;
    void set(size_t, size_t, double);
};
Matrix triplet_matrix(const std::vector<Tree> &, size_t, const Config &, bool allow_gpu = true);
std::vector<Id> sorted_rows(const Matrix &, const Config &);
void complete_tree(Tree &, const std::vector<Id> &, size_t, bool);
Tree upgma(Matrix, const Config &);
void add_s4_candidates(CandidateSpace &, const std::vector<Tree> &, const Registry &,
                       const Config &, Scope = Scope::LINEAGE);
std::vector<Tree> generate_s6(const std::vector<Tree> &, const Registry &, const Config &);
bool cuda_available(std::string *reason = nullptr);
std::vector<Split> cuda_cross_candidates(const CandidateSpace &, const Config &);
std::vector<Score> cuda_weights(const CandidateSpace &, const Observations &,
                                const std::vector<Id> &, const Config &);
bool cuda_triplet_matrix(const std::vector<Tree> &, size_t, const Config &, Matrix &);
struct Focus {
    size_t n = 0;
    std::vector<Id> selected;
    std::vector<uint64_t> pair, sister;
    uint64_t support(Id, Id, Id, Id focus) const;
};
class CudaFocus {
  public:
    virtual ~CudaFocus() = default;
    virtual Focus focus(Id, const std::vector<Id> &) = 0;
};
std::unique_ptr<CudaFocus> make_cuda_focus(const Observations &, size_t, const Config &);
Focus cpu_focus(const Observations &, size_t, Id, const std::vector<Id> &);
} // namespace stelar
