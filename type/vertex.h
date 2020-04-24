#pragma once

#include <list>

#include "common/register.h"
#include "matrix/dense.h"
#include "geometry/jet.h"

namespace hitnlls {

class VertexBase;

template <typename TupleType, int i>
struct VertexTypeAt {
    using VertexType = typename VertexTypeAt<typename TupleType::BaseType, i - 1>::VertexType;
};
template <typename TupleType>
struct VertexTypeAt<TupleType, 0> {
    using VertexType = typename TupleType::VertexType;
};
template <typename TupleType, int i>
struct VertexAt {
    static inline typename VertexTypeAt<TupleType, i>::VertexType *Vertex(TupleType &tuple) {
        return VertexAt<typename TupleType::BaseType, i - 1>::Vertex(tuple);
    }
};
template <typename TupleType>
struct VertexAt<TupleType, 0> {
    static inline typename TupleType::VertexType *Vertex(TupleType &tuple) {
        return tuple.vptr_;
    }
};
template <typename ...V> class VertexTuple;
template <> class VertexTuple<> { template <typename TupleType, int i> friend struct VertexAt; };
template <typename V> 
class VertexTuple<V> {
public:
    static constexpr int VertexNums = 1;
    using ThisType = VertexTuple<V>;
    using VertexType = V;

    template <typename TupleType, int i>
    friend struct VertexAt;

    VertexTuple() { vptr_ = nullptr; }
    VertexTuple(V *vptr) { vptr_ = vptr; }

    void SetVertices(V *vptr) { vptr_ = vptr; }

    template <int i>
    typename VertexTypeAt<ThisType, i>::VertexType *GetVertex() {
        return VertexAt<ThisType, i>::Vertex(*this);
    }
    const VertexBase *GetVertex(int i) const {
        if (i == 0) {
            return vptr_;
        } else {
            return nullptr;
        }
    }
    VertexBase *GetVertex(int i) {
        if (i == 0) {
            return vptr_;
        } else {
            return nullptr;
        }
    }
    int GetVertexId(int i) const {
        if (i == 0) {
            return vptr_->GetId();
        } else {
            -1;
        }
    }

    static constexpr int TupleSize() { return 1; }
    template <int i>
    static constexpr int PerturbationDim() {
        return (i == 0) ? V::PerturbationDim : 0;
    }
    static constexpr int PerturbationDim(int i) {
        return (i == 0) ? V::PerturbationDim : 0;
    }
    template <int i>
    static constexpr int PerturbationOffset() {
        return (i == 0) ? 0 : V::PerturbationDim;
    }
    static constexpr int PerturbationOffset(int i) {
        return (i == 0) ? 0 : V::PerturbationDim;
    }
    static constexpr int TotalPerturbationDim = PerturbationOffset<VertexNums>();

private:
    V *vptr_;
};
template <typename V, typename ...VRest>
class VertexTuple<V, VRest...> : private VertexTuple<VRest...> {
public:
    static constexpr int VertexNums = sizeof...(VRest) + 1;
    using ThisType = VertexTuple<V, VRest...>;
    using BaseType = VertexTuple<VRest...>;
    using VertexType = V;

    template <typename TupleType, int i>
    friend struct VertexAt;

    VertexTuple() { vptr_ = nullptr; }
    VertexTuple(V *vptr, VRest*... vrest_ptrs) : BaseType(vrest_ptrs...) { vptr_ = vptr; }

    void SetVertices(V *vptr, VRest*... vrest_ptrs) { BaseType::SetVertices(vrest_ptrs...); vptr_ = vptr; }

    template <int i>
    typename VertexTypeAt<ThisType, i>::VertexType *GetVertex() {
        return VertexAt<ThisType, i>::Vertex(*this);
    }
    const VertexBase *GetVertex(int i) const {
        if (i == 0) {
            return vptr_;
        } else {
            return BaseType::GetVertex(i - 1);
        }
    }
    VertexBase *GetVertex(int i) {
        if (i == 0) {
            return vptr_;
        } else {
            return BaseType::GetVertex(i - 1);
        }
    }
    int GetVertexId(int i) const {
        if (i == 0) {
            return vptr_->GetId();
        } else {
            return BaseType::GetVertexId(i - 1);
        }
    }

    static constexpr int TupleSize() { return sizeof...(VRest) + 1; }
    template <int i>
    static constexpr int PerturbationDim() {
        return (i == 0) ? V::PerturbationDim : BaseType::template PerturbationDim<i - 1>();
    }
    static constexpr int PerturbationDim(int i) {
        return (i == 0) ? V::PerturbationDim : BaseType::PerturbationDim(i - 1);
    }
    template <int i>
    static constexpr int PerturbationOffset() {
        return (i == 0) ? 0 : BaseType::template PerturbationOffset<i - 1>() + V::PerturbationDim;
    }
    static constexpr int PerturbationOffset(int i) {
        return (i == 0) ? 0 : BaseType::PerturbationOffset(i - 1) + V::PerturbationDim;
    }
    static constexpr int TotalPerturbationDim = PerturbationOffset<VertexNums>();

private:
    V *vptr_;
};

class VertexBase {
public:
    enum Status { Active = 0, NonActive, Fixed };

    explicit VertexBase() { id_ = 0; status_ = Active; marginalized_ = false; }
    virtual ~VertexBase() = default;

    void SetId(int id) { id_ = id; }
    void SetStatus(Status status) { status_ = status; }
    void SetMarginalized() { marginalized_ = true; }

    int GetId() { return id_; }
    Status GetStatus() { return status_; }
    bool GetMarginalized() { return marginalized_; }

    virtual int GetPerturbationDim() const = 0;

    virtual void SetZero() = 0;
    virtual void ApplyPerturbation(const double *inc) = 0;
    virtual void Push() = 0;
    virtual void Pop() = 0;

    VertexBase(const VertexBase &) = delete;
    VertexBase &operator=(const VertexBase &) = delete;
    VertexBase &operator=(const VertexBase &&) = delete;

protected:
    int id_;
    Status status_;
    bool marginalized_;
};

template <int PDIM, typename EstimateType>
class Vertex : public VertexBase {
public:
    static constexpr int PerturbationDim = PDIM;
    using PerturbationVector = matrix::Matrix<double, PDIM, 1>;

    const EstimateType &GetEstimate() const { return estimate_; }
    virtual int GetPerturbationDim() const override final { return PDIM; }
    virtual void ApplyPerturbation(const double *inc) override {
        matrix::Map<const PerturbationVector> pert(inc);
        ApplyPerturbation(pert);
    }
    virtual void Push() override { stack_.push_back(estimate_); }
    virtual void Pop() override { SetEstimate(stack_.back()); stack_.pop_back(); }
    virtual void SetEstimate(const EstimateType &est) { estimate_ = est; }

    virtual void SetZero() = 0;
    virtual void ApplyPerturbation(const PerturbationVector &pert) = 0;

protected:
    EstimateType estimate_;
    std::list<EstimateType> stack_;
};

template <int PDIM, template <typename Scalar> class GenericEstimateType>
class VertexAudoDiff : public Vertex<PDIM, GenericEstimateType<double>> {
public:
    using ThisType = VertexAudoDiff<PDIM, GenericEstimateType>;
    using BaseType = Vertex<PDIM, GenericEstimateType<double>>;
    using PerturbationVector = matrix::Matrix<double, PDIM, 1>;
    using EstimateType = GenericEstimateType<double>;
    using PerturbationJetVector = matrix::Matrix<geometry::Jetd, PDIM, 1>;
    using EstimateTypeJet = GenericEstimateType<geometry::Jetd>;

    const EstimateTypeJet &GetJetEstimate() const { return estimate_jet_; }
    virtual void Pop() override { SetEstimate(BaseType::stack_.back()); BaseType::stack_.pop_back(); }
    virtual void ApplyPerturbation(const PerturbationVector &pert) override {
        PerturbationJetVector pert_jet;
        geometry::ConvertMatrix(pert_jet, pert);
        ApplyPerturbation(pert_jet);
        geometry::ConvertMatrix(BaseType::estimate_, estimate_jet_);
    }
    virtual void SetEstimate(const EstimateType &est) override {
        BaseType::estimate_ = est;
        geometry::ConvertMatrix(estimate_jet_, BaseType::estimate_);
    }

    virtual void SetZero() = 0;
    virtual void ApplyPerturbation(const PerturbationJetVector &pert) = 0;

protected:
    EstimateTypeJet estimate_jet_;
};

HITNLLS_REGISTER_REGISTER(VertexBase)
#define HITNLLS_REGISTER_VERTEX(name) \
    HITNLLS_REGISTER_CLASS(VertexBase, name)

} // namespace hitnlls