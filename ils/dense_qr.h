#pragma once

#include <unordered_map>

#include "solver.h"

namespace nlls {
namespace internal {
class DenseQr : public Solver {
public:
    explicit DenseQr(Strategy s) : s_(s) {}
protected:
    virtual void BuildStructure() override;
    virtual void FillStructure() override;
    virtual float InitLevenMarqParam() override;
    virtual VectorXf ComputeGaussNewtonPoint(float l = 0.0) override;
    virtual VectorXf ComputeCauchyPoint() override;
    virtual void Update(const VectorXf &inc) override;
    virtual float SurrogateReduction(const VectorXf &inc) override;

    Strategy s_;
    MatrixXf A_;
    VectorXf b_;
    std::unordered_map<VertexBase *, int> v2i_;
    std::unordered_map<FactorBase *, int> f2i_;
};
} // namespace internal
} // namespace nlls