#pragma once

#include "matrix/dense.h"
#include "type/vertex.h"

namespace hitnlls {

template <typename T> using Vector3 = matrix::Matrix<T, 3, 1>;

struct VertexPt3 : public Vertex<3, Vector3<double>> {
    virtual void SetZero() override;
    virtual void ApplyPerturbation(const matrix::Vector3d &pert) override;
};

struct VertexPt3Jet : public VertexAudoDiff<3, Vector3> {
    using BaseType = VertexAudoDiff<3, Vector3>;
    using PerturbationJetVector = typename BaseType::PerturbationJetVector;

    virtual void SetZero() {
        matrix::Vector3d v;
        v.SetZero();
        SetEstimate(v);
    }
    virtual void ApplyPerturbation(const PerturbationJetVector &pert)  {
        estimate_jet_ += pert;
    }
};

} // namespace hitnlls