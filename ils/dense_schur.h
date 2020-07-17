#pragma once

#include <unordered_map>

#include "solver.h"

namespace nlls {
namespace internal {

class DenseSchur : public Solver {
public:
    explicit DenseSchur(Strategy s) : s_(s) {}

protected:
    virtual void BuildStructure() override;
    virtual void FillStructure() override;
    virtual float InitLevenMarqParam() override;
    virtual VectorXf ComputeGaussNewtonPoint(float l = 0.0) override;
    virtual VectorXf ComputeCauchyPoint() override;
    virtual void Update(const VectorXf &inc) override;
    virtual float SurrogateReduction(const VectorXf &inc) override;

    float WeightedNormSquared(const VectorXf &inc);

    Strategy s_;
    std::unordered_map<VertexBase *, int> vr2i_;
    std::unordered_map<VertexBase *, int> vs2i_;
    std::unordered_map<VertexBase *, MatrixXf> vs2c_;
    std::unordered_map<VertexBase *, std::unordered_map<VertexBase *, MatrixXf>> vs2vr2c_;
    MatrixXf Ar_;
    VectorXf br_;
    VectorXf bs_;
};

} // namespace internal
} // namespace nlls