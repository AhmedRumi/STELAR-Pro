#include "stelar/core.hpp"
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#ifndef STELAR_SOURCE_ROOT
#define STELAR_SOURCE_ROOT "."
#endif
namespace fs = std::filesystem;
namespace {
using namespace stelar;
constexpr const char *version = "STELAR-Pro C++17 0.1.0";
fs::path repository_root() {
    const char *override_root = std::getenv("STELAR_PRO_HOME");
    return override_root && *override_root ? fs::path(override_root) : fs::path(STELAR_SOURCE_ROOT);
}
struct TempFile {
    std::string path;
    explicit TempFile(std::string prefix = "/tmp/stelar-cpp-") {
        prefix += "XXXXXX";
        std::vector<char> b(prefix.begin(), prefix.end());
        b.push_back(0);
        int fd = mkstemp(b.data());
        if (fd < 0)
            throw std::runtime_error("Cannot create temporary file");
        close(fd);
        path = b.data();
    }
    ~TempFile() {
        std::error_code ec;
        fs::remove(path, ec);
    }
    TempFile(const TempFile &) = delete;
};
void write_file(const std::string &path, const std::string &text) {
    if (path.empty()) {
        std::cout << text;
        return;
    }
    auto parent = fs::absolute(path).parent_path();
    fs::create_directories(parent);
    TempFile tmp((parent / ".stelar-cpp-").string());
    {
        std::ofstream f(tmp.path);
        f << text;
        if (!f)
            throw std::runtime_error("Cannot write " + path);
    }
    fs::rename(tmp.path, path);
}
std::string read_file(const std::string &p) {
    std::ifstream f(p);
    if (!f)
        throw std::runtime_error("Cannot read " + p);
    return {std::istreambuf_iterator<char>(f), {}};
}
int run_process(const std::vector<std::string> &args, bool capture_output = true) {
    TempFile capture;
    pid_t pid = fork();
    if (pid < 0)
        throw std::runtime_error("Cannot start subprocess");
    if (pid == 0) {
        int fd = capture_output ? open(capture.path.c_str(), O_WRONLY | O_TRUNC) : -1;
        if (fd >= 0) {
            dup2(fd, 1);
            dup2(fd, 2);
            close(fd);
        }
        std::vector<char *> a;
        for (auto &s : args)
            a.push_back(const_cast<char *>(s.c_str()));
        a.push_back(nullptr);
        execvp(a[0], a.data());
        _exit(127);
    }
    int status;
    while (waitpid(pid, &status, 0) < 0)
        if (errno != EINTR)
            throw std::runtime_error("waitpid failed");
    if (!WIFEXITED(status) || WEXITSTATUS(status))
        throw std::runtime_error("Backend failed: " + args.front() + "\n" +
                                 read_file(capture.path));
    return WEXITSTATUS(status);
}
void tag_input(const Config &cfg, const std::string &input, const std::string &out) {
    std::vector<std::string> args{cfg.astral_executable, "-T"};
    if (!cfg.gene_map.empty()) {
        args.push_back("-a");
        args.push_back(cfg.gene_map);
    }
    args.insert(args.end(), {"-o", out, input});
    run_process(args);
    if (fs::file_size(out) == 0)
        throw std::runtime_error("Rooting/tagging backend produced no trees");
}
class Tee : public std::streambuf {
    std::streambuf *old;
    std::ofstream file;

  public:
    explicit Tee(const std::string &p) : old(std::cerr.rdbuf()), file(p) {
        if (!file)
            throw std::runtime_error("Cannot open log: " + p);
        std::cerr.rdbuf(this);
    }
    ~Tee() {
        std::cerr.rdbuf(old);
    }
    int overflow(int c) override {
        if (c != EOF) {
            old->sputc(char(c));
            file.put(char(c));
        }
        return c;
    }
    int sync() override {
        old->pubsync();
        file.flush();
        return 0;
    }
};
void help() {
    std::cout
        << version
        << "\nUsage: stelar-pro-cpp -i genes.trees [-o species.tree] [options]\n"
           "  --search-space s1|s2|s3   Java S1, S4, S6 respectively (default s1)\n"
           "  --cpu | --gpu | --auto | --gpu-strict\n"
           "  -t, --threads N          CPU worker limit\n"
           "  --cross-tree-splits | --no-cross-tree-splits\n"
           "  --search-mode local|full  Independent cross-tree override\n"
           "  --exhaustive             s2: add S1 and globally complete DISCO\n"
           "  --s3-rounds N, --s6-rounds N   s3 guide rounds (default 4)\n"
           "  --s3-seed N, --s6-seed N       Java-compatible search seed\n"
           "  --max-search-passes 1    One outer pass\n"
           "  --taxa-file FILE         Infer/score induced species universe\n"
           "  -c, --score-species-tree FILE  Score a fixed species tree\n"
           "  -T, --tag-only           Root/tag and write -o\n"
           "  --extract-taxa --taxa-set union|intersection\n"
           "  --astral-pro-executable FILE  Root/tag backend override\n"
           "  --gene-species-map FILE  Mapping for root/tag backend\n"
           "  --no-prune-search-space | --prune-search-space\n"
           "  --large-n-score-type int128|double\n"
           "  --seeds N               Compatibility hash setting (sets also checked exactly)\n"
           "  --no-gpu-batch | --gpu-batch-size N | --gpu-batches N\n"
           "  --gpu-vram-occupancy-factor F | --gpu-vram-control-factor F\n"
           "  --gpu-treewalk-vram-cap-mb N | --gpu-sim-vram-cap-mb N\n"
           "  --gpu-progress-interval F | --gpu-dist-tile-size N\n"
           "  --gpu-dp-state-space-construction-output-cap SIZE\n"
           "  --gpu-dp-state-space-progress-time-interval F\n"
           "  --gpu-dp-state-space-progress-max-steps N\n"
           "  --dump-clusters FILE | --dump-completed-gene-trees FILE\n"
           "  --verify-parse | --verify-hash | --verify-clusters | --verify-partitions\n"
           "  --verify-dp | --verify-weights | --verify-similarity-matrix | --verify-upgma\n"
           "  --verify-distance-matrix  CPU diagnostic, independent of inference\n"
           "  --verify-greedy-consensus  Legacy verifier via unchanged Java sources (JDK needed)\n"
           "  --log-file FILE | -vv | -vvv | -q | --quiet\n"
           "  --rooted | --keep-polytomy-during-inference\n"
           "  --completion-method similarity|distance (legacy selector)\n"
           "  --stepb-restriction dlogd|n | --stepb-fast-restriction\n"
           "  --stepb-quadratic-nn-balls | --stepb-random-leftover-resolution\n"
           "  --stepb-process-large-polytomies\n"
           "  --no-build | --no-notify | --xms SIZE | --xmx SIZE (launcher compatibility)\n"
           "  --diagnose | --version | --help\n"
           "  --tagged-input          Diagnostic: input is already rooted and event-tagged\n"
           "s3 uses the s2 pipeline at <=100 species. Consensus/autocomplete and anchored\n"
           "unrooted configurations retain Java's unsupported-configuration errors.\n";
}
long long integer(const std::string &s) {
    size_t n;
    auto v = std::stoll(s, &n);
    if (n != s.size())
        throw std::invalid_argument("Invalid integer: " + s);
    return v;
}
double real(const std::string &s) {
    size_t n;
    auto v = std::stod(s, &n);
    if (n != s.size() || !std::isfinite(v))
        throw std::invalid_argument("Invalid number: " + s);
    return v;
}
Config parse(int argc, char **argv) {
    Config c;
    c.threads = std::max(1U, std::thread::hardware_concurrency());
    std::set<std::string> values = {"--search-space",
                                    "--search-mode",
                                    "--s3-rounds",
                                    "--s6-rounds",
                                    "--s3-seed",
                                    "--s6-seed",
                                    "--max-search-passes",
                                    "--threads",
                                    "--num-threads",
                                    "-t",
                                    "--seeds",
                                    "-m",
                                    "--input",
                                    "-i",
                                    "--output",
                                    "-o",
                                    "--log-file",
                                    "--score-species-tree",
                                    "--species-tree",
                                    "--score",
                                    "-c",
                                    "--taxa-file",
                                    "--species-list",
                                    "--species-list-file",
                                    "--taxa-set",
                                    "--taxa-operation",
                                    "--astral-pro-executable",
                                    "--gene-species-map",
                                    "--large-n-score-type",
                                    "--large-score-type",
                                    "--anchor-taxon",
                                    "--gpu-batch-size",
                                    "--gpu-batches",
                                    "--gpu-vram-occupancy-factor",
                                    "--gpu-vram-control-factor",
                                    "--gpu-treewalk-vram-cap-mb",
                                    "--gpu-progress-interval",
                                    "--gpu-dist-tile-size",
                                    "--gpu-sim-vram-cap-mb",
                                    "--gpu-dp-state-space-construction-output-cap",
                                    "--gpu-dp-state-space-progress-time-interval",
                                    "--gpu-dp-state-space-progress-max-steps",
                                    "--dump-clusters",
                                    "--dump-completed-gene-trees",
                                    "--completion-method",
                                    "--stepb-restriction",
                                    "--xms",
                                    "--Xms",
                                    "--xmx",
                                    "--Xmx"};
    std::set<std::string> switches = {"--cpu",
                                      "--gpu",
                                      "--auto",
                                      "--gpu-strict",
                                      "--cross-tree-splits",
                                      "--no-cross-tree-splits",
                                      "--exhaustive",
                                      "--no-prune-search-space",
                                      "--no-prune-unreachable",
                                      "--prune-search-space",
                                      "--prune-unreachable",
                                      "--rooted",
                                      "--unrooted",
                                      "--anchor-outgroup",
                                      "--anchor",
                                      "--no-anchor-outgroup",
                                      "--no-anchor",
                                      "--keep-polytomy-during-inference",
                                      "--no-gpu-batch",
                                      "--verify-parse",
                                      "--verify-hash",
                                      "--verify-clusters",
                                      "--verify-partitions",
                                      "--verify-dp",
                                      "--verify-weights",
                                      "--verify-distance-matrix",
                                      "--verify-similarity-matrix",
                                      "--verify-upgma",
                                      "--verify-greedy-consensus",
                                      "--consensus-experimental",
                                      "--stepb-fast-restriction",
                                      "--stepb-quadratic-nn-balls",
                                      "--stepb-random-leftover-resolution",
                                      "--stepb-process-large-polytomies",
                                      "--resolve-input-gene-tree-polytomies",
                                      "--autocomplete-incomplete-gene-trees",
                                      "--extract-taxa",
                                      "-T",
                                      "--tag-only",
                                      "-v",
                                      "-vv",
                                      "-vvv",
                                      "-q",
                                      "--quiet",
                                      "--no-build",
                                      "--no-notify",
                                      "-nn",
                                      "--diagnose",
                                      "--tagged-input"};
    int position = 0;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i], v = "1";
        auto eq = k.find('=');
        if (eq != std::string::npos) {
            v = k.substr(eq + 1);
            k.resize(eq);
        }
        if (k == "--help" || k == "-h") {
            help();
            std::exit(0);
        }
        if (k == "--version" || k == "-v") {
            std::cout << version << '\n';
            std::exit(0);
        }
        if (k.empty() || k[0] != '-') {
            if (position++ == 0)
                c.input = k;
            else if (position == 2)
                c.output = k;
            else
                throw std::invalid_argument("Unexpected positional argument: " + k);
            continue;
        }
        if (values.count(k)) {
            if (eq == std::string::npos) {
                if (i + 1 == argc)
                    throw std::invalid_argument(k + " requires a value");
                v = argv[++i];
            }
            if (v.empty())
                throw std::invalid_argument(k + " requires a value");
        } else if (!switches.count(k))
            throw std::invalid_argument("Unknown option: " + k);
        c.options[k] = v;
        if (k == "--input" || k == "-i")
            c.input = v;
        else if (k == "--output" || k == "-o")
            c.output = v;
        else if (k == "--log-file")
            c.log_file = v;
        else if (k == "--search-space") {
            // Keep preset reset order compatible with Java CliPresets; cross
            // tree overrides remain independent of preset position.
            for (auto reset :
                 {"--autocomplete-incomplete-gene-trees", "--consensus-experimental",
                  "--resolve-input-gene-tree-polytomies", "--stepb-quadratic-nn-balls",
                  "--stepb-random-leftover-resolution", "--stepb-process-large-polytomies"})
                c.options.erase(reset);
            c.completion_method = "similarity";
            if (v.size() > 1 && (v[0] == 's' || v[0] == 'S'))
                v.erase(0, 1);
            c.mode = integer(v);
            if (c.mode < 1 || c.mode > 3)
                throw std::invalid_argument(
                    "Search spaces are s1 (Java S1), s2 (Java S4), s3 (Java S6)");
        } else if (k == "--cpu")
            c.compute = Config::CPU;
        else if (k == "--gpu")
            c.compute = Config::GPU;
        else if (k == "--auto")
            c.compute = Config::AUTO;
        else if (k == "--gpu-strict") {
            c.compute = Config::GPU;
            c.gpu_strict = true;
        } else if (k == "--threads" || k == "--num-threads" || k == "-t")
            c.threads = integer(v);
        else if (k == "--seeds" || k == "-m")
            c.seeds = integer(v);
        else if (k == "--s3-rounds" || k == "--s6-rounds")
            c.rounds = integer(v);
        else if (k == "--s3-seed" || k == "--s6-seed")
            c.seed = integer(v);
        else if (k == "--cross-tree-splits")
            c.cross_tree = true;
        else if (k == "--no-cross-tree-splits")
            c.cross_tree = false;
        else if (k == "--search-mode") {
            if (v != "local" && v != "full")
                throw std::invalid_argument("--search-mode expects local|full");
            c.cross_tree = v == "full";
        } else if (k == "--exhaustive")
            c.exhaustive = true;
        else if (k == "--no-prune-search-space" || k == "--no-prune-unreachable")
            c.prune = false;
        else if (k == "--prune-search-space" || k == "--prune-unreachable")
            c.prune = true;
        else if (k == "--score-species-tree" || k == "--species-tree" || k == "--score" ||
                 k == "-c") {
            c.score_only = true;
            c.score_tree = v;
        } else if (k == "--taxa-file" || k == "--species-list" || k == "--species-list-file")
            c.taxa_file = v;
        else if (k == "--gene-species-map")
            c.gene_map = v;
        else if (k == "--astral-pro-executable")
            c.astral_executable = v;
        else if (k == "--large-n-score-type" || k == "--large-score-type") {
            if (v == "exact" || v == "i128")
                v = "int128";
            if (v == "fp64" || v == "float")
                v = "double";
            if (v != "int128" && v != "double")
                throw std::invalid_argument("Expected int128|double");
            c.score_type = v;
        } else if (k == "--no-gpu-batch")
            c.no_gpu_batch = true;
        else if (k == "--gpu-batch-size")
            c.gpu_batch_size = integer(v);
        else if (k == "--gpu-batches")
            c.gpu_batches = integer(v);
        else if (k == "--gpu-vram-occupancy-factor")
            c.gpu_occupancy = real(v);
        else if (k == "--gpu-vram-control-factor")
            c.gpu_control = real(v);
        else if (k == "--gpu-treewalk-vram-cap-mb")
            c.gpu_treewalk_cap_mb = integer(v);
        else if (k == "--gpu-sim-vram-cap-mb")
            c.gpu_sim_cap_mb = integer(v);
        else if (k == "--gpu-dist-tile-size")
            c.gpu_tile = integer(v);
        else if (k == "--gpu-progress-interval")
            c.gpu_progress = real(v);
        else if (k == "--gpu-dp-state-space-progress-time-interval")
            c.gpu_dp_progress = real(v);
        else if (k == "--gpu-dp-state-space-progress-max-steps")
            c.gpu_dp_steps = integer(v);
        else if (k == "--gpu-dp-state-space-construction-output-cap") {
            double factor = 1;
            if (!v.empty() && std::isalpha(v.back())) {
                char u = std::tolower(v.back());
                factor = u == 'k' ? 1e3 : u == 'm' ? 1e6 : u == 'g' ? 1e9 : 0;
                v.pop_back();
            }
            double bytes = real(v) * factor;
            if (!factor || bytes < 12 || bytes / 12 > INT32_MAX)
                throw std::invalid_argument("Invalid GPU output cap");
            c.gpu_dp_cap = bytes / 12;
        } else if (k == "--dump-clusters")
            c.dump_clusters = v;
        else if (k == "--dump-completed-gene-trees")
            c.dump_completed = v;
        else if (k == "--completion-method") {
            if (v != "similarity" && v != "distance")
                throw std::invalid_argument("Completion method expects similarity|distance");
            c.completion_method = v;
        } else if (k == "--stepb-restriction") {
            if (v != "dlogd" && v != "n" && v != "fast" && v != "full")
                throw std::invalid_argument("--stepb-restriction expects dlogd|n");
        } else if (k == "--taxa-set" || k == "--taxa-operation") {
            if (v != "union" && v != "intersection")
                throw std::invalid_argument("Taxa operation expects union|intersection");
        } else if (k == "--keep-polytomy-during-inference")
            c.keep_polytomy = true;
        else if (k == "--no-anchor" || k == "--no-anchor-outgroup") {
            c.options.erase("--anchor");
            c.options.erase("--anchor-outgroup");
        } else if (k == "--rooted")
            c.options.erase("--unrooted");
        else if (k == "-q" || k == "--quiet") {
            c.quiet = true;
            c.verbosity = 0;
        } else if (k == "-vv")
            c.verbosity = 2;
        else if (k == "-vvv")
            c.verbosity = 3;
    }
    if (c.threads < 1 || c.seeds < 1 || c.rounds < 1 || c.gpu_batch_size < 0 || c.gpu_batches < 0 ||
        c.gpu_sim_cap_mb < 1 || c.gpu_treewalk_cap_mb < 1 || c.gpu_tile < 0 ||
        c.gpu_occupancy <= 0 || c.gpu_occupancy > 1 || c.gpu_control < 0 || c.gpu_progress < 0 ||
        c.gpu_dp_steps < 1 || c.gpu_dp_progress < 0)
        throw std::invalid_argument("Invalid positive count or GPU memory limit");
    if (c.input.empty() && !c.option("--diagnose"))
        throw std::invalid_argument("-i/--input is required");
    if (c.option("--max-search-passes") && integer(c.options.at("--max-search-passes")) != 1)
        throw std::invalid_argument("Retained modes support one outer pass");
    if (c.exhaustive && c.mode != 2)
        throw std::invalid_argument("--exhaustive applies to s2 (Java S4)");
    if (c.mode == 1 && c.cross_tree)
        throw std::invalid_argument("s1 uses local transitions (same restriction as Java S1)");
    if (c.option("--unrooted") || c.option("--anchor") || c.option("--anchor-outgroup"))
        throw std::invalid_argument(
            "STELAR-Pro uses a rooted objective; anchoring/unrooted inference is unsupported");
    for (auto k : {"--consensus-experimental", "--resolve-input-gene-tree-polytomies",
                   "--autocomplete-incomplete-gene-trees"})
        if (c.option(k))
            throw std::invalid_argument(
                std::string(k) + ": unsupported search configuration, matching production Java");
    if (c.option("--extract-taxa") && (c.score_only || !c.taxa_file.empty()))
        throw std::invalid_argument(
            "--extract-taxa cannot be combined with scoring or a taxa file");
    if ((c.option("-T") || c.option("--tag-only")) &&
        (c.score_only || c.option("--extract-taxa") || c.output.empty()))
        throw std::invalid_argument(
            "Tag-only requires output and cannot combine other utility modes");
    if (!c.gene_map.empty() && !c.option("-T") && !c.option("--tag-only") &&
        !c.option("--extract-taxa") && !c.option("--diagnose"))
        throw std::invalid_argument(
            "Gene-copy mapping ingestion is unsupported for inference/scoring, as in Java. "
            "Use species-mapped input. --gene-species-map remains available for --tag-only.");
    if (c.astral_executable.empty()) {
        const char *p = std::getenv("STELAR_PRO_EXECUTABLE");
        c.astral_executable = p ? p : (repository_root() / "ASTER-Linux/bin/astral-pro3").string();
    }
    return c;
}
void validate_paths(const Config &);
// Consensus/Step A/Step B are outside all three retained search algorithms.
// Preserve their old diagnostic flag using the unchanged reference verifier,
// instead of silently running a different consensus algorithm in production.
int legacy_consensus_verifier(const Config &cfg, int argc, char **argv) {
    validate_paths(cfg);
    std::vector<char> path;
    const std::string prefix = "/tmp/stelar-cpp-consensus-verifier-XXXXXX";
    path.assign(prefix.begin(), prefix.end());
    path.push_back(0);
    if (!mkdtemp(path.data()))
        throw std::runtime_error("Cannot create legacy verifier classes directory");
    struct Cleanup {
        std::string path;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(path, ec); // Only this mkdtemp-created directory.
        }
    } cleanup{path.data()};
    log(cfg, "Legacy consensus diagnostic: compiling/running the unchanged Java reference. "
             "Native s1/s2/s3 inference does not use Java or consensus.");
    const fs::path source_root = repository_root();
    run_process({"javac", "-sourcepath", (source_root / "src").string(), "-d", cleanup.path,
                 (source_root / "src/stelarx/Main.java").string()});
    std::vector<std::string> command{"java",
                                     "-Dstelarpro.home=" + source_root.string(),
                                     "-Djava.library.path=" + (source_root / "native").string(),
                                     "-cp",
                                     cleanup.path,
                                     "stelarx.Main"};
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i], key = arg.substr(0, arg.find('='));
        if (key == "--no-build" || key == "--no-notify" || key == "-nn")
            continue;
        if (key == "--xms" || key == "--Xms" || key == "--xmx" || key == "--Xmx") {
            if (arg.find('=') == std::string::npos)
                ++i;
            continue;
        }
        if (key == "--tagged-input")
            throw std::invalid_argument(
                "Legacy Java consensus verifier has no --tagged-input mode");
        if (key == "--search-space") {
            std::string value =
                arg.find('=') == std::string::npos ? argv[++i] : arg.substr(arg.find('=') + 1);
            int mode =
                integer(value.front() == 's' || value.front() == 'S' ? value.substr(1) : value);
            command.insert(command.end(), {"--search-space", mode == 1   ? "S1"
                                                             : mode == 2 ? "S4"
                                                                         : "S6"});
        } else {
            if (key == "--s3-rounds")
                arg.replace(0, key.size(), "--s6-rounds");
            if (key == "--s3-seed")
                arg.replace(0, key.size(), "--s6-seed");
            command.push_back(arg);
        }
    }
    return run_process(command, false);
}
void validate_paths(const Config &c) {
    std::vector<std::string> reads = {c.input, c.score_tree, c.taxa_file, c.gene_map,
                                      c.astral_executable};
    std::vector<std::string> writes = {c.output, c.log_file, c.dump_clusters, c.dump_completed};
    if (!c.output.empty())
        writes.push_back(c.output + ".stats.json");
    if (c.option("--verify-weights"))
        writes.push_back(c.output + ".weights.tsv");
    auto same = [](const std::string &a, const std::string &b) {
        return fs::weakly_canonical(fs::absolute(a)) == fs::weakly_canonical(fs::absolute(b)) ||
               (fs::exists(a) && fs::exists(b) && fs::equivalent(a, b));
    };
    for (size_t i = 0; i < writes.size(); ++i)
        if (!writes[i].empty()) {
            for (auto &r : reads)
                if (!r.empty() && same(writes[i], r))
                    throw std::invalid_argument("Output aliases an input: " + writes[i]);
            for (size_t j = 0; j < i; ++j)
                if (!writes[j].empty() && same(writes[i], writes[j]))
                    throw std::invalid_argument("Output paths must differ");
        }
}
void filter_taxa(std::vector<Tree> &trees, Registry &r, const Config &cfg) {
    if (cfg.taxa_file.empty())
        return;
    std::istringstream in(read_file(cfg.taxa_file));
    std::string name;
    Registry selected;
    std::vector<Id> ids;
    Registry species_registry;
    if (cfg.score_only)
        read_trees(cfg.score_tree, species_registry);
    while (std::getline(in, name)) {
        size_t a = name.find_first_not_of(" \t\r");
        if (a == std::string::npos || name[a] == '#')
            continue;
        size_t b = name.find_last_not_of(" \t\r");
        name = name.substr(a, b - a + 1);
        Id i = r.find(name);
        if (i != none && selected.find(name) == none &&
            (!cfg.score_only || species_registry.find(name) != none)) {
            selected.add(name);
            ids.push_back(i);
        }
    }
    if (ids.size() < 3)
        throw std::invalid_argument(
            "Fewer than three listed taxa occur in the effective tree universe");
    trees = restrict_genes(trees, ids, r.size());
    trees.erase(std::remove_if(trees.begin(), trees.end(),
                               [&](auto &t) {
                                   return cfg.score_only ? cardinality(t.taxa()) < 3
                                                         : t.leaves.size() < 2;
                               }),
                trees.end());
    r = std::move(selected);
}
void dump_sets(const CandidateSpace &c, const Registry &r, const std::string &path) {
    std::ostringstream out;
    out << "cluster_id\tsize\ttaxa\n";
    for (Id i = 0; i < c.eligible.size(); ++i)
        if (c.eligible[i]) {
            out << i << '\t' << c.sets.count(i) << '\t';
            bool first = true;
            for (Id t : members(c.sets.bits(i))) {
                if (!first)
                    out << ',';
                first = false;
                out << r.names[t];
            }
            out << '\n';
        }
    write_file(path, out.str());
}
// Legacy S1 matrix diagnostics do not participate in inference. In s2/s3,
// similarity diagnostics instead expose their observed DISCO triplet evidence.
Matrix diagnostic_matrix(const std::vector<Tree> &trees, size_t n, bool distance) {
    Matrix result(n), denominator(n);
    auto add = [&](Id a, Id b, double num, double den) {
        if (a == b)
            return;
        result.set(a, b, result.get(a, b) + num);
        denominator.set(a, b, denominator.get(a, b) + den);
    };
    auto choose2 = [](size_t x) { return x > 1 ? double(x) * double(x - 1) / 2 : 0.; };
    for (const auto &t : trees) {
        std::vector<Id> depth(t.nodes.size()), leaf_depth(n);
        std::vector<Id> stack{t.root};
        while (!stack.empty()) {
            Id u = stack.back();
            stack.pop_back();
            auto ch = t.children(u);
            if (ch.empty())
                leaf_depth[t.nodes[u].taxon] = depth[u];
            for (auto it = ch.rbegin(); it != ch.rend(); ++it) {
                depth[*it] = depth[u] + 1;
                stack.push_back(*it);
            }
        }
        size_t k = t.leaves.size();
        if (!distance && k < 4)
            continue;
        for (Id u : t.postorder()) {
            auto ch = t.children(u);
            if (ch.empty())
                continue;
            if (ch.size() != 2)
                throw std::invalid_argument("Legacy matrix diagnostic requires binary input");
            auto a = t.nodes[ch[0]], b = t.nodes[ch[1]], v = t.nodes[u];
            size_t outside = k - (v.hi - v.lo);
            for (size_t i = a.lo; i < a.hi; ++i)
                for (size_t j = b.lo; j < b.hi; ++j) {
                    Id x = t.leaves[i], y = t.leaves[j];
                    add(x, y,
                        distance ? double(leaf_depth[x]) + leaf_depth[y] - 2. * depth[u]
                                 : choose2(outside),
                        distance ? 1. : 0.);
                }
            if (!distance) {
                for (size_t j = 0; j < k; ++j)
                    if (j < v.lo || j >= v.hi) {
                        for (size_t i = a.lo; i < a.hi; ++i)
                            add(t.leaves[i], t.leaves[j], choose2(b.hi - b.lo), 0.);
                        for (size_t i = b.lo; i < b.hi; ++i)
                            add(t.leaves[i], t.leaves[j], choose2(a.hi - a.lo), 0.);
                    }
            }
        }
        if (!distance)
            for (size_t i = 0; i < k; ++i)
                for (size_t j = i + 1; j < k; ++j)
                    add(t.leaves[i], t.leaves[j], 0., choose2(k - 2));
    }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            result.set(i, j,
                       denominator.get(i, j) ? result.get(i, j) / denominator.get(i, j)
                       : distance            ? std::numeric_limits<double>::infinity()
                                             : 0.);
    return result;
}
bool matrix_diagnostic(const std::vector<Tree> &genes, const Registry &reg, const Config &cfg) {
    bool distance = cfg.option("--verify-distance-matrix");
    if (!distance && !cfg.option("--verify-similarity-matrix") && !cfg.option("--verify-upgma"))
        return false;
    Matrix matrix;
    if (distance || cfg.mode == 1) {
        log(cfg, "Legacy matrix diagnostic: CPU; not an inference stage");
        matrix = diagnostic_matrix(genes, reg.size(), distance);
    } else {
        auto disco = decompose(genes, reg.size(),
                               cfg.exhaustive  ? Scope::GLOBAL
                               : cfg.mode == 3 ? Scope::DUPLICATION
                                               : Scope::LINEAGE,
                               cfg);
        log(cfg,
            "Triplet diagnostic over observed DISCO components (s3 has subset-local matrices)");
        matrix = triplet_matrix(disco.trees, reg.size(), cfg);
    }
    std::ostringstream out;
    out << (distance                       ? "DISTANCE_MATRIX"
            : cfg.option("--verify-upgma") ? "UPGMA_BIPARTITIONS"
                                           : "SIMILARITY_MATRIX")
        << "\nn=" << reg.size() << "\ntaxa=";
    for (size_t i = 0; i < reg.size(); ++i)
        out << (i ? "," : "") << reg.names[i];
    out << '\n';
    if (!distance && cfg.option("--verify-upgma")) {
        auto guide = upgma(std::move(matrix), cfg);
        for (Id u : guide.postorder())
            if (u != guide.root && guide.nodes[u].child != none) {
                std::vector<std::string> names;
                for (size_t i = guide.nodes[u].lo; i < guide.nodes[u].hi; ++i)
                    names.push_back(reg.names[guide.leaves[i]]);
                std::sort(names.begin(), names.end());
                out << "bipartition=";
                for (size_t i = 0; i < names.size(); ++i)
                    out << (i ? "," : "") << names[i];
                out << '\n';
            }
    } else {
        out << std::fixed << std::setprecision(distance ? 6 : 8);
        for (size_t i = 0; i < reg.size(); ++i) {
            out << (distance ? "row" : "sim_row") << i << '=';
            for (size_t j = 0; j < reg.size(); ++j)
                out << (j ? "," : "") << (i == j ? distance ? 0. : 1. : matrix.get(i, j));
            out << '\n';
        }
    }
    write_file(cfg.output, out.str());
    return true;
}
Score fixed_score(const Tree &t, const Observations &o, const Registry &r, const Config &cfg) {
    CandidateSpace c(r.size());
    c.add_tree(t, false);
    c.finalize();
    auto ids = c.reachable(false);
    Score s = 0;
    std::vector<Score> w;
    if (cfg.compute == Config::GPU)
        w = cuda_weights(c, o, ids, cfg);
    else {
        w.resize(ids.size());
        parallel_for(ids.size(), cfg.threads, [&](size_t j) {
            auto &sp = c.splits[ids[j]];
            w[j] = split_weight(c.sets, sp.left, sp.right, o);
        });
    }
    for (auto v : w)
        s += v;
    return s;
}
int execute(Config cfg) {
    validate_paths(cfg);
    std::unique_ptr<Tee> tee;
    if (!cfg.log_file.empty()) {
        fs::create_directories(fs::absolute(cfg.log_file).parent_path());
        tee = std::make_unique<Tee>(cfg.log_file);
    }
    std::string gpu;
    bool available = cuda_available(&gpu);
    if (cfg.compute != Config::CPU) {
        if (!available && cfg.gpu_strict)
            throw std::runtime_error("--gpu-strict: " + gpu);
        cfg.compute = available ? Config::GPU : Config::CPU;
    }
    log(cfg, std::string(version) + "; " + (cfg.compute == Config::GPU ? "CUDA: " + gpu : "CPU") +
                 "; " + std::to_string(cfg.threads) + " workers");
    if (cfg.score_type == "double")
        log(cfg, "--large-n-score-type double accepted for CLI compatibility; native weights "
                 "remain exact (64-bit when safe, otherwise 128-bit), without FP64 rounding.");
    if (cfg.option("--diagnose")) {
        std::cout << "CUDA: " << gpu << "\nRoot/tag executable: " << cfg.astral_executable
                  << "\nNode bytes: " << sizeof(Node) << "\nJNI: none\n";
        return 0;
    }
    auto begin = std::chrono::steady_clock::now();
    if (cfg.option("--extract-taxa")) {
        Registry r;
        auto trees = read_trees(cfg.input, r);
        Bits selected = trees[0].taxa();
        std::string op = "union";
        if (cfg.option("--taxa-set"))
            op = cfg.options["--taxa-set"];
        if (cfg.option("--taxa-operation"))
            op = cfg.options["--taxa-operation"];
        if (op != "union" && op != "intersection")
            throw std::invalid_argument("Taxa operation expects union|intersection");
        for (size_t i = 1; i < trees.size(); ++i)
            selected = op == "union" ? set_union(selected, trees[i].taxa())
                                     : set_intersection(selected, trees[i].taxa());
        std::vector<std::string> names;
        for (Id i : members(selected))
            names.push_back(r.names[i]);
        std::sort(names.begin(), names.end());
        std::string out;
        for (auto &n : names)
            out += n + "\n";
        write_file(cfg.output, out);
        return 0;
    }
    if (cfg.option("-T") || cfg.option("--tag-only")) {
        TempFile tagged;
        tag_input(cfg, cfg.input, tagged.path);
        write_file(cfg.output, read_file(tagged.path));
        return 0;
    }
    Registry input_registry;
    auto input_trees = read_trees(cfg.input, input_registry);
    bool polytomy = false;
    for (auto &t : input_trees)
        if (!t.binary())
            polytomy = true;
    input_trees.clear();
    input_trees.shrink_to_fit();
    TempFile resolved, tagged;
    std::string preprocessing = cfg.input;
    if (!cfg.score_only && !cfg.option("--tagged-input")) {
        const char *py = std::getenv("STELAR_PRO_PYTHON");
        run_process({py ? py : "python3",
                     (repository_root() / "scripts/arb_resolve_polytomies.py").string(), cfg.input,
                     resolved.path});
        preprocessing = resolved.path;
    }
    if (!cfg.option("--tagged-input")) {
        log(cfg, "Preprocessing: ASTRAL-Pro rooting/tagging");
        tag_input(cfg, preprocessing, tagged.path);
        preprocessing = tagged.path;
    }
    Registry reg;
    auto genes = read_trees(preprocessing, reg);
    filter_taxa(genes, reg, cfg);
    polytomy = false;
    for (auto &t : genes)
        if (!t.binary())
            polytomy = true;
    if (reg.size() < 1)
        throw std::invalid_argument("No species remain");
    log(cfg, "Parsed " + std::to_string(genes.size()) + " gene families and " +
                 std::to_string(reg.size()) + " species");
    if (cfg.option("--verify-parse")) {
        std::string out;
        for (auto &t : genes)
            out += newick(t, reg, true) + "\n";
        write_file(cfg.output, out);
        return 0;
    }
    if (cfg.score_only) {
        Registry species_reg = cfg.taxa_file.empty() ? reg : Registry{};
        auto list = read_trees(cfg.score_tree, species_reg, !cfg.taxa_file.empty());
        if (list.size() != 1)
            throw std::invalid_argument("Expected one species tree");
        Tree t = std::move(list.front());
        if (!cfg.taxa_file.empty()) {
            std::vector<int> map(species_reg.size(), -1);
            for (Id i = 0; i < species_reg.size(); ++i) {
                Id j = reg.find(species_reg.names[i]);
                if (j != none)
                    map[i] = j;
            }
            t = restrict_tree(t, map);
        }
        if (t.leaves.size() != reg.size() || cardinality(t.taxa()) != reg.size())
            throw std::invalid_argument(
                "Species tree must contain every effective taxon exactly once");
        auto observations = build_observations(genes);
        auto score = fixed_score(t, observations, reg, cfg);
        std::string text = "Final triplet score = " + score_string(score) + "\n";
        log(cfg, text);
        write_file(cfg.output, text);
        return 0;
    }
    if (cfg.mode == 3 && reg.size() <= 100) {
        log(cfg, "s3: <=100 species; using s2 (Java S4)");
        cfg.mode = 2;
    }
    if (matrix_diagnostic(genes, reg, cfg))
        return 0;
    auto prep = std::chrono::steady_clock::now();
    CandidateSpace candidates(reg.size());
    if (cfg.mode == 1) {
        for (auto &t : genes)
            candidates.add_tree(t, true);
    } else if (cfg.mode == 2)
        add_s4_candidates(candidates, genes, reg, cfg);
    else {
        auto guides = generate_s6(genes, reg, cfg);
        for (auto &t : guides)
            candidates.add_tree(t, false);
    }
    double prep_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - prep).count();
    log(cfg, "Candidate preparation: " + std::to_string(prep_seconds) + " s");
    if (cfg.cross_tree) {
        size_t before = candidates.splits.size();
        candidates.add_cross(cfg);
        log(cfg, "Cross-tree splits: +" + std::to_string(candidates.splits.size() - before));
    }
    if (!cfg.dump_clusters.empty())
        dump_sets(candidates, reg, cfg.dump_clusters);
    if (cfg.option("--verify-hash") || cfg.option("--verify-clusters"))
        dump_sets(candidates, reg, cfg.output.empty() ? "cpp-clusters.tsv" : cfg.output);
    auto observations = build_observations(genes);
    log(cfg,
        "Original unique rooted partitions: " + std::to_string(observations.partitions.size()));
    if (cfg.option("--verify-partitions")) {
        std::ostringstream f;
        for (auto &p : observations.partitions) {
            f << p.frequency;
            for (size_t j = 0; j < p.children; ++j) {
                f << '\t';
                for (Id tax :
                     members(observations.sets.bits(observations.child_sets[p.offset + j])))
                    f << reg.names[tax] << ',';
            }
            f << '\n';
        }
        write_file(cfg.output, f.str());
        return 0;
    }
    if (cfg.option("--verify-dp")) {
        std::ostringstream f;
        for (auto &s : candidates.splits) {
            for (Id id : {s.parent, s.left, s.right}) {
                for (Id tax : members(candidates.sets.bits(id)))
                    f << reg.names[tax] << ',';
                f << '\t';
            }
            f << '\n';
        }
        write_file(cfg.output, f.str());
        return 0;
    }
    if (cfg.option("--verify-hash") || cfg.option("--verify-clusters"))
        return 0;
    auto result = infer(candidates, observations, reg, cfg);
    std::string tree = newick(result.tree, reg);
    write_file(cfg.output, tree + "\n");
    Score final_score = result.score;
    if (polytomy && !cfg.keep_polytomy && !cfg.option("--tagged-input")) {
        log(cfg, "Final scoring against unresolved input polytomies");
        TempFile raw_tagged;
        tag_input(cfg, cfg.input, raw_tagged.path);
        Registry rawreg;
        auto rawgenes = read_trees(raw_tagged.path, rawreg);
        std::vector<int> map(rawreg.size(), -1);
        for (Id i = 0; i < rawreg.size(); ++i) {
            Id j = reg.find(rawreg.names[i]);
            if (j != none)
                map[i] = j;
        }
        for (auto &t : rawgenes)
            t = restrict_tree(t, map);
        auto rawobs = build_observations(rawgenes);
        final_score = fixed_score(result.tree, rawobs, reg, cfg);
    }
    log(cfg, "Inference DP: optimization-objective triplet score = " + score_string(result.score));
    log(cfg, "Final triplet score = " + score_string(final_score));
    double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    log(cfg, "Total time: " + std::to_string(seconds) +
                 " s; peak RSS: " + std::to_string(usage.ru_maxrss) + " KiB");
    if (!cfg.output.empty()) {
        std::ostringstream s;
        s << "{\n  \"mode\": \"s" << cfg.mode << "\",\n  \"taxa\": " << reg.size()
          << ",\n  \"families\": " << genes.size() << ",\n  \"clusters\": " << result.clusters
          << ",\n  \"splits\": " << result.splits << ",\n  \"scored_splits\": " << result.scored
          << ",\n  \"partitions\": " << observations.partitions.size()
          << ",\n  \"optimization_score\": \"" << score_string(result.score)
          << "\",\n  \"final_score\": \"" << score_string(final_score)
          << "\",\n  \"seconds\": " << std::setprecision(12) << seconds
          << ",\n  \"candidate_seconds\": " << prep_seconds
          << ",\n  \"peak_rss_kib\": " << usage.ru_maxrss << "\n}\n";
        write_file(cfg.output + ".stats.json", s.str());
    }
    return 0;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 1) {
            help();
            return 2;
        }
        Config cfg = parse(argc, argv);
        if (cfg.option("--verify-greedy-consensus") && !cfg.score_only &&
            !cfg.option("--extract-taxa") && !cfg.option("--verify-parse") && !cfg.option("-T") &&
            !cfg.option("--tag-only"))
            return legacy_consensus_verifier(cfg, argc, argv);
        return execute(std::move(cfg));
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 2;
    }
}
