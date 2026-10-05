// Host-side S5 bridge: the shared A-Pro placement/NNI code invokes the CUDA
// per-locus quartet counter through the same ordered scorer interface as CPU.
#define STELAR_S5_GPU_BUILD
#include "stelar_pro_s5_greedy.cpp"
#include "stelar_pro_s5_cuda_api.h"

#include <deque>

namespace {

class CudaQuartetBackend final : public QuartetScoreBackend {
    void* context = nullptr;
    std::vector<S5CudaOperation> pending;
    std::deque<score_t> ready;
    int scoreCount = 0;

    void flush() {
        if (pending.empty()) throw std::logic_error("S5 CUDA scorer queue is empty");
        std::vector<long long> sums(scoreCount);
        s5CudaExecute(context, static_cast<int>(pending.size()), pending.data(),
                      scoreCount, sums.data());
        for (const auto& op : pending) {
            ready.push_back(op.color == -2 ? sums[op.scoreSlot] : 0);
        }
        pending.clear();
        scoreCount = 0;
    }

    score_t pop() {
        if (ready.empty()) flush();
        score_t value = ready.front();
        ready.pop_front();
        return value;
    }

public:
    CudaQuartetBackend(const TripartitionInitializer& input, int taxa) {
        const int genes = static_cast<int>(input.nodes.size());
        std::vector<S5CudaNode> nodes;
        std::vector<int> offsets;
        std::vector<int> leaves;
        offsets.reserve(static_cast<size_t>(genes) * taxa + 1);
        for (int gene = 0; gene < genes; ++gene) {
            const int base = static_cast<int>(nodes.size());
            for (const auto& source : input.nodes[gene]) {
                S5CudaNode node{};
                node.up = source.up < 0 ? -1 : base + source.up;
                node.small = source.small < 0 ? -1 : base + source.small;
                node.large = source.large < 0 ? -1 : base + source.large;
                node.dup = source.dup ? 1 : 0;
                nodes.push_back(node);
            }
            for (int taxon = 0; taxon < taxa; ++taxon) {
                offsets.push_back(static_cast<int>(leaves.size()));
                if (taxon < static_cast<int>(input.leafParent[gene].size())) {
                    for (int local : input.leafParent[gene][taxon]) {
                        leaves.push_back(base + local);
                    }
                }
            }
        }
        offsets.push_back(static_cast<int>(leaves.size()));
        context = s5CudaCreate(genes, taxa, static_cast<int>(nodes.size()),
                               nodes.data(), offsets.data(),
                               static_cast<int>(leaves.size()), leaves.data());
    }

    ~CudaQuartetBackend() override { s5CudaDestroy(context); }

    void updateSet(int color, int taxon) override {
        pending.push_back({color, taxon, -1});
    }
    void updateGet() override { (void)pop(); }
    void scoreSet() override {
        pending.push_back({-2, -1, scoreCount++});
    }
    score_t scoreGet() override { return pop(); }
};

}  // namespace

std::unique_ptr<QuartetScoreBackend> makeS5GpuBackend(
    const TripartitionInitializer& genes, int taxa) {
    return std::unique_ptr<QuartetScoreBackend>(new CudaQuartetBackend(genes, taxa));
}
