#pragma once

#include <vector>
#include <unordered_map>

#include "solver.h"

namespace nlls {
namespace internal {

class SparseCholesky : public Solver {
public:
    explicit SparseCholesky(Strategy s) : s_(s) {}

protected:
    virtual void BuildStructure() override;
    virtual void FillStructure() override;
    virtual float InitLevenMarqParam() override;
    virtual VectorXf ComputeGaussNewtonPoint(float l = 0.0) override;
    virtual VectorXf ComputeCauchyPoint() override;
    virtual void Update(const VectorXf &inc) override;
    virtual float SurrogateReduction(const VectorXf &inc) override;

    void ComputeOrdering();
    float WeightedNormSquared(const VectorXf &inc);

    Strategy s_;
    std::vector<VertexBase *> orders_;
    std::unordered_map<VertexBase *, int> v2i_;
    std::unordered_map<VertexBase *, std::unordered_map<VertexBase *, MatrixXf>> vc2vr2c_;
    VectorXf b_;
};

} // namespace internal
} // namespace nlls