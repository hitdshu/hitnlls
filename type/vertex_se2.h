#pragma once

#include "../ils/vertex.h"
#include "../matrix/dense.h"

namespace nlls {

template <typename T> using Matrix33 = Matrix<T, 3, 3>;

struct VertexSE2 : public Vertex<3, Matrix33f> {
    virtual void Reset() override;
    virtual void ApplyPerturbation(const Vector3f &pert) override;
};

struct VertexSE2Jet : public VertexAutoDiff<3, Matrix33> {
    virtual void Reset() override;
    virtual void ApplyPerturbation(const Vector3j &pert) override;
};

} // namespace nlls