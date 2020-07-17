#pragma once

#include "../ils/vertex.h"
#include "../matrix/dense.h"

namespace nlls {

template <typename T> using Matrix44 = Matrix<T, 4, 4>;

struct VertexSE3 : public Vertex<6, Matrix44f> {
    virtual void Reset() override;
    virtual void ApplyPerturbation(const Vector6f &pert) override;
};

struct VertexSE3Jet : public VertexAutoDiff<6, Matrix44> {
    virtual void Reset() override;
    virtual void ApplyPerturbation(const Vector6j &pert) override;
};

} // namespace nlls