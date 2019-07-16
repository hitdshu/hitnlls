#pragma once

#include "graph/base_graph.h"

namespace hitnlls {
namespace solver {

class BaseSolver {
public:
    BaseSolver() { graph_ = nullptr; steps_ = 10; error_reduction_ = 1e-4; early_stop_ = true; }
    virtual ~BaseSolver() = default;

    virtual void SetGraph(::hitnlls::graph::BaseGraph *graph) final { graph_ = graph; }
    virtual void SetOptimizationSteps(int steps) final { steps_ = steps; }
    virtual void SetErrorReduction(float error_reduction) final { error_reduction_ = error_reduction; }
    virtual void SetEarlyStop(bool early_stop) final { early_stop_ = early_stop; }

    virtual void Optimize() = 0;

    BaseSolver(const BaseSolver &) = delete;
    BaseSolver &operator=(const BaseSolver &) = delete;

protected:
    ::hitnlls::graph::BaseGraph *graph_;
    int steps_;
    float error_reduction_;
    bool early_stop_;
};

HITNLLS_REGISTER_REGISTERER(BaseSolver);
#define HITNLLS_REGISTER_SOLVER(name) \
    HITNLLS_REGISTER_CLASS(BaseSolver, name)

} // namespace solver
} // namespace hitnlls