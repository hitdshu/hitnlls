#pragma once

#include <unordered_set>

#include "../common/macros.h"

#include "factor.h"
#include "vertex.h"

namespace nlls {
namespace internal {
struct LevenMarqparam { 
    float u; 
    float v; 
    explicit LevenMarqparam() : u(-1), v(2) {}
};

struct TrustRegionParam {
    float r;
    float v;
    float l;
    float u;
    explicit TrustRegionParam() : r(-1), v(2), l(0.25), u(0.75) {}
};
} // namespace internal

enum Strategy { TRUST_REGION, LEVEN_MARQ, GAUSS_NEWTON };

enum Method { DENSE_QR, DENSE_CHOLESKY, DENSE_SCHUR, SPARSE_CHOLESKY };

struct SolverOption {
    int num_iters;
    int num_inner_iters;
    float term_ratio;
    float term_error;
    Strategy sttg;
    Method method;
    explicit SolverOption() : num_iters(10), num_inner_iters(10), term_ratio(1e-5), term_error(1e-10), sttg(TRUST_REGION), method(DENSE_QR) {}
};

class Problem;

namespace internal {

class Solver {
public:
    NLLS_NONCOPYABLE(Solver)
    explicit Solver() = default;
    virtual ~Solver() = default;

    float ComputeChi() const;
    std::unordered_set<FactorBase *> &Factors();
    const std::unordered_set<FactorBase *> &Factors() const;
    std::unordered_set<VertexBase *> &Vertices();
    const std::unordered_set<VertexBase *> &Vertices() const;
    std::unordered_set<VertexBase *> &SchurVertices();
    const std::unordered_set<VertexBase *> &SchurVertices() const;

    void SetProblem(Problem *p) { problem_ = p; }
    void Solve();

    static Solver *CreateSolver(Problem *p, Strategy s, Method m);

protected:
    void PushVertices();
    void PopVertices();

    float SolveOneStep();
    float SolveOneStepGaussNewton();
    float SolveOneStepLevenMarq();
    float SolveOneStepTrustRegion();

    virtual void BuildStructure() = 0;
    virtual void FillStructure() = 0;
    virtual float InitLevenMarqParam() = 0;
    virtual VectorXf ComputeGaussNewtonPoint(float l = 0.0) = 0;
    virtual VectorXf ComputeCauchyPoint() = 0;
    virtual void Update(const VectorXf &inc) = 0;
    virtual float SurrogateReduction(const VectorXf &inc) = 0;

    Problem *problem_;
    
    internal::LevenMarqparam lmp_;
    internal::TrustRegionParam trp_;
};

} // namespace internal
} // namespace nlls
