#include "stelar/core.hpp"

namespace stelar {
bool cuda_available(std::string *reason) {
    if (reason)
        *reason = "this executable was built without CUDA support";
    return false;
}
std::vector<Split> cuda_cross_candidates(const CandidateSpace &, const Config &) {
    throw std::runtime_error("CUDA cross-tree construction is unavailable in this executable");
}
std::vector<Score> cuda_weights(const CandidateSpace &, const Observations &,
                                const std::vector<Id> &, const Config &) {
    throw std::runtime_error("CUDA weight calculation is unavailable in this executable");
}
bool cuda_triplet_matrix(const std::vector<Tree> &, size_t, const Config &, Matrix &) {
    return false;
}
std::unique_ptr<CudaFocus> make_cuda_focus(const Observations &, size_t, const Config &) {
    return nullptr;
}
} // namespace stelar
