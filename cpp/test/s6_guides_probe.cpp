#include "../src/s6.cpp"
#include <fstream>
#include <numeric>
int main(int argc, char **argv) {
    using namespace stelar;
    if (argc != 3)
        return 2;
    Registry registry;
    auto genes = read_trees(argv[1], registry);
    Config cfg;
    cfg.compute = Config::CPU;
    cfg.threads = 1;
    cfg.quiet = true;
    auto disco = decompose(genes, registry.size(), Scope::DUPLICATION, cfg);
    std::vector<Id> order(registry.size());
    std::iota(order.begin(), order.end(), 0);
    JavaRandom random(cfg.seed);
    random.shuffle(order);
    size_t n = registry.size(),
           m = std::max<size_t>(4, std::min(n, size_t(std::sqrt(double(n)) * std::log(double(n)) /
                                                      std::log(2.0) / 4.0)));
    order.resize(m);
    auto guides = backbone_guides(disco.trees, order, n, cfg);
    Registry subset;
    for (Id taxon : order)
        subset.add(registry.names[taxon]);
    std::ofstream out(argv[2]);
    for (auto &tree : guides)
        out << newick(tree, subset) << '\n';
    CandidateSpace candidates(m);
    for (auto &tree : guides)
        candidates.add_tree(tree);
    candidates.finalize();
    auto ids = candidates.reachable(true);
    std::vector<std::string> lines;
    auto names = [&](Id id) {
        std::vector<std::string> list;
        for (Id t : members(candidates.sets.bits(id)))
            list.push_back(subset.names[t]);
        std::sort(list.begin(), list.end());
        std::string s;
        for (auto &name : list) {
            if (!s.empty())
                s += ",";
            s += name;
        }
        return s;
    };
    for (Id id : ids) {
        auto &s = candidates.splits[id];
        auto a = names(s.left), b = names(s.right);
        lines.push_back(a < b ? a + " | " + b : b + " | " + a);
    }
    std::sort(lines.begin(), lines.end());
    std::ofstream splits(std::string(argv[2]) + ".splits");
    for (auto &line : lines)
        splits << line << '\n';
}
