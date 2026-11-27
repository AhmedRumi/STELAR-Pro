#include "stelar/core.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace stelar;

static void ensure(bool ok, const char *what) {
    if (!ok)
        throw std::runtime_error(what);
}

static void oracle() {
    Registry registry;
    auto a = parse_newick("((A,B),(C,D));", registry, true, false);
    auto b = parse_newick("((A,C),B);", registry, true, false);
    Config cfg;
    cfg.compute = Config::CPU;
    cfg.threads = 4;
    cfg.quiet = true;
    const auto matrix = triplet_matrix({a, b}, 4, cfg, false);
    const Id A = registry.find("A"), B = registry.find("B"), C = registry.find("C"),
             D = registry.find("D");
    ensure(matrix.get(A, B) == 2.0 / 3.0, "AB triplet matrix oracle");
    ensure(matrix.get(A, C) == 1.0 / 3.0, "AC triplet matrix oracle");
    ensure(matrix.get(C, D) == 1, "CD triplet matrix oracle");
    ensure(matrix.get(A, D) == 0, "AD triplet matrix oracle");
    Config single = cfg;
    single.threads = 1;
    ensure(triplet_matrix({a, b}, 4, single, false).values == matrix.values,
           "matrix thread reproducibility");
    Matrix zero(6);
    auto rows = sorted_rows(zero, cfg);
    for (Id taxon = 0; taxon < 6; ++taxon) {
        ensure(rows[taxon * 6] == taxon, "self distance rank");
        Id rank = 1;
        for (Id other = 6; other-- > 0;)
            if (other != taxon)
                ensure(rows[taxon * 6 + rank++] == other, "descending ID distance tie");
    }
    auto tree = upgma(Matrix(6), cfg);
    Registry numbers;
    for (int i = 0; i < 6; ++i)
        numbers.add(std::to_string(i));
    ensure(newick(tree, numbers) == "(((((0,1),5),4),3),2);", "UPGMA cached-neighbor equal ties");
    ensure(newick(upgma(Matrix(6), single), numbers) == newick(tree, numbers),
           "UPGMA thread reproducibility");
    std::cout << "S4 hand-computed triplet, anchor and UPGMA tie oracles passed\n";
}

// This dump is compared by a Python harness with JavaS4Reference, preserving
// numeric IDs, child order, scope masks and exact double values.
static void reference(const std::string &input, const std::string &output, int threads) {
    Registry registry;
    const auto genes = read_trees(input, registry, true, false);
    Config cfg;
    cfg.compute = Config::CPU;
    cfg.threads = threads;
    cfg.quiet = true;
    std::ofstream out(output);
    if (!out)
        throw std::runtime_error("reference dump output");
    out << "N " << registry.size() << '\n';
    for (Id x = 0; x < registry.size(); ++x)
        out << "T " << x << ' ' << registry.names[x] << '\n';
    const Scope modes[] = {Scope::GLOBAL, Scope::LINEAGE, Scope::DUPLICATION};
    const char *names[] = {"GLOBAL", "LINEAGE", "DUPLICATION"};
    for (size_t mode = 0; mode < 3; ++mode) {
        auto d = decompose(genes, registry.size(), modes[mode], cfg);
        out << "S " << names[mode] << ' ' << d.cuts << ' ' << d.discarded << ' ' << d.trees.size()
            << '\n';
        if (d.trees.empty())
            continue;
        for (const auto &t : d.trees) {
            out << "R " << newick(t, registry) << " M";
            for (Id x : members(t.missing))
                out << ' ' << x;
            out << '\n';
        }
        const auto matrix = triplet_matrix(d.trees, registry.size(), cfg, false);
        const auto rows = sorted_rows(matrix, cfg);
        for (Id a = 0; a < registry.size(); ++a) {
            for (Id b = a + 1; b < registry.size(); ++b)
                out << "P " << a << ' ' << b << ' ' << std::hexfloat << matrix.get(a, b)
                    << std::defaultfloat << '\n';
            out << "A " << a;
            for (Id rank = 0; rank < registry.size(); ++rank)
                out << ' ' << rows[a * registry.size() + rank];
            out << '\n';
        }
        out << "U " << newick(upgma(matrix, cfg), registry) << '\n';
        for (auto &t : d.trees) {
            complete_tree(t, rows, registry.size(), modes[mode] != Scope::GLOBAL);
            out << "C " << newick(t, registry) << '\n';
        }
    }
}

int main(int argc, char **argv) {
    try {
        oracle();
        if (argc >= 3)
            reference(argv[1], argv[2], argc >= 4 ? std::stoi(argv[3]) : 4);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "S4 test failed: " << e.what() << '\n';
        return 1;
    }
}
