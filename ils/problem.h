#pragma once

#include <unordered_set>
#include <vector>
#include <string>

#include "../common/macros.h"

#include "vertex.h"
#include "factor.h"
#include "solver.h"

namespace nlls {

struct IterStat {
    int iter;
    float chi;
    float inc;
};

class Problem {
public:
    NLLS_NONCOPYABLE(Problem)
    explicit Problem() = default;
    ~Problem();

    void SetOption(const SolverOption &option) { option_ = option; }
    const SolverOption &GetOption() const { return option_; }
    const IterStat &GetInitStat() const { return *stat_.begin(); }
    const IterStat &GetFinalStat() const { return stat_.back(); }
    const std::vector<IterStat> &GetAllStat() const { return stat_; }
    std::string BriefReport() const;
    std::string Report() const;

    void AddFactor(FactorBase *f) { factors_.insert(f); }
    void AddSchurVertex(VertexBase *v) { if (v->GetStatus() == VertexBase::Active) { schur_vertices_.insert(v); } }
    float Solve();

private:
    void BuildStructure();

    std::unordered_set<FactorBase *> factors_;
    std::unordered_set<VertexBase *> vertices_;
    std::unordered_set<VertexBase *> schur_vertices_;
    std::vector<IterStat> stat_;
    SolverOption option_;
    friend class internal::Solver;
};

} // namespace nlls