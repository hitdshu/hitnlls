#pragma once

#include "nlls/graph.h"

namespace hitnlls {

class SolverBase {
public:
    SolverBase() { graph_ = nullptr; steps_ = 10; error_reduction_ = 1e-4; early_stop_ = true; }
    virtual ~SolverBase() = default;

    void SetGraph(GraphBase *graph) { graph_ = graph; }
    void SetOptimizationSteps(int steps) { steps_ = steps; }
    void SetErrorReduction(double error_reduction) { error_reduction_ = error_reduction; }
    void SetEarlyStop(bool early_stop) { early_stop_ = early_stop; }

    virtual void Optimize() = 0;

    SolverBase(const SolverBase &) = delete;
    SolverBase &operator=(const SolverBase &) = delete;
    SolverBase &operator=(SolverBase &&) = delete;

protected:
    GraphBase *graph_;
    int steps_;
    double error_reduction_;
    bool early_stop_;
};

HITNLLS_REGISTER_REGISTER(SolverBase)
#define HITNLLS_REGISTER_SOLVER(name) \
    HITNLLS_REGISTER_CLASS(SolverBase, name)

} // namespace hitnlls