#pragma once

#include "../ils/vertex.h"

namespace nlls {

template <int N>
class VertexRn : public Vertex<N> {
public:
    using BaseType = Vertex<N>;
    using PerturbationVector = typename BaseType::PerturbationVector;
    virtual void Reset() override {
        BaseType::estimate_.SetZero();
    }
    virtual void ApplyPerturbation(const PerturbationVector &pert) override {
        BaseType::estimate_ += pert;
    }
};

template <typename T> using Vector1 = Matrix<T, 1, 1>;
class VertexR1Jet : public VertexAutoDiff<1, Vector1> {
public:
    using BaseType = VertexAutoDiff<1, Vector1>;
    using PerturbationVectorJet = typename BaseType::PerturbationVectorJet;
    virtual void Reset() override {
        BaseType::estimate_jet_.SetZero();
    }
    virtual void ApplyPerturbation(const PerturbationVectorJet &pert) override {
        BaseType::estimate_jet_ += pert;
    }
};

template <typename T> using Vector2 = Matrix<T, 2, 1>;
class VertexR2Jet : public VertexAutoDiff<2, Vector2> {
public:
    using BaseType = VertexAutoDiff<2, Vector2>;
    using PerturbationVectorJet = typename BaseType::PerturbationVectorJet;
    virtual void Reset() override {
        BaseType::estimate_jet_.SetZero();
    }
    virtual void ApplyPerturbation(const PerturbationVectorJet &pert) override {
        BaseType::estimate_jet_ += pert;
    }
};

template <typename T> using Vector3 = Matrix<T, 3, 1>;
class VertexR3Jet : public VertexAutoDiff<3, Vector3> {
public:
    using BaseType = VertexAutoDiff<3, Vector3>;
    using PerturbationVectorJet = typename BaseType::PerturbationVectorJet;
    virtual void Reset() override {
        BaseType::estimate_jet_.SetZero();
    }
    virtual void ApplyPerturbation(const PerturbationVectorJet &pert) override {
        BaseType::estimate_jet_ += pert;
    }
};

template <typename T> using Vector4 = Matrix<T, 4, 1>;
class VertexR4Jet : public VertexAutoDiff<4, Vector4> {
public:
    using BaseType = VertexAutoDiff<4, Vector4>;
    using PerturbationVectorJet = typename BaseType::PerturbationVectorJet;
    virtual void Reset() override {
        BaseType::estimate_jet_.SetZero();
    }
    virtual void ApplyPerturbation(const PerturbationVectorJet &pert) override {
        BaseType::estimate_jet_ += pert;
    }
};

template <typename T> using Vector5 = Matrix<T, 5, 1>;
class VertexR5Jet : public VertexAutoDiff<5, Vector5> {
public:
    using BaseType = VertexAutoDiff<5, Vector5>;
    using PerturbationVectorJet = typename BaseType::PerturbationVectorJet;
    virtual void Reset() override {
        BaseType::estimate_jet_.SetZero();
    }
    virtual void ApplyPerturbation(const PerturbationVectorJet &pert) override {
        BaseType::estimate_jet_ += pert;
    }
};

template <typename T> using Vector6 = Matrix<T, 6, 1>;
class VertexR6Jet : public VertexAutoDiff<6, Vector6> {
public:
    using BaseType = VertexAutoDiff<6, Vector6>;
    using PerturbationVectorJet = typename BaseType::PerturbationVectorJet;
    virtual void Reset() override {
        BaseType::estimate_jet_.SetZero();
    }
    virtual void ApplyPerturbation(const PerturbationVectorJet &pert) override {
        BaseType::estimate_jet_ += pert;
    }
};

} // namespace nlls