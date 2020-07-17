#pragma once

#include <list>

#include "../common/macros.h"
#include "../common/traits.h"
#include "../common/register.h"
#include "../matrix/dense.h"

namespace nlls {

class VertexBase {
public:
    NLLS_NONCOPYABLE(VertexBase)
    enum Status { Active, NonActive, Fixed };
    explicit VertexBase() { id_ = 0; status_ = Active; }
    virtual ~VertexBase() = default;
    NLLS_INLINE void SetId(int id) { id_ = id; }
    NLLS_INLINE int GetId() { return id_; }
    NLLS_INLINE void SetStatus(Status status) { status_ = status; }
    NLLS_INLINE Status GetStatus() { return status_; }
    virtual int Dof() const = 0;
    virtual void Reset() = 0;
    virtual void Update(const float *inc) = 0;
    virtual void Push() = 0;
    virtual void Pop() = 0;
private:
    int id_;
    Status status_;
};
template <int N, typename T = Matrix<float, N, 1>>
class Vertex : public VertexBase {
public:
    static constexpr int DOF = N;
    using EstimateType = T;
    using PerturbationVector = Matrix<float, DOF, 1>;
    virtual int Dof() const override final { return DOF; }
    virtual void Update(const float *inc) override { ApplyPerturbation(Map<const PerturbationVector>(inc)); }
    virtual void Push() override { stack_.push_back(estimate_); }
    virtual void Pop() override { SetEstimate(stack_.back()); stack_.pop_back(); }
    NLLS_INLINE const EstimateType &GetEstimate() const { return estimate_; }
    virtual void SetEstimate(const EstimateType &est) { estimate_ = est; }
    virtual void Reset() = 0;
    virtual void ApplyPerturbation(const PerturbationVector &pert) = 0;
protected:
    EstimateType estimate_;
    std::list<EstimateType> stack_;
};
template <int N, template <typename T> class GenericEstimateType>
class VertexAutoDiff : public Vertex<N, GenericEstimateType<float>> {
public:
    using ThisType = VertexAutoDiff;
    using BaseType = Vertex<N, GenericEstimateType<float>>;
    using PerturbationVector = Matrix<float, N, 1>;
    using EstimateType = GenericEstimateType<float>;
    using PerturbationVectorJet = Matrix<Jetf, N, 1>;
    using EstimateTypeJet = GenericEstimateType<Jetf>;
    virtual void SetEstimate(const EstimateType &est) override { BaseType::SetEstimate(est); ConvertMatrix(estimate_jet_, est); }
    virtual void ApplyPerturbation(const PerturbationVector &pert) override final {
        PerturbationVectorJet pert_jet;
        ConvertMatrix(pert_jet, pert);
        ApplyPerturbation(pert_jet);
        ConvertMatrix(BaseType::estimate_, estimate_jet_);
    }
    NLLS_INLINE const EstimateTypeJet &GetEstimateJet() const { return estimate_jet_; }
    virtual void Reset() = 0;
    virtual void ApplyPerturbation(const PerturbationVectorJet &pert_jet) = 0;
protected:
    EstimateTypeJet estimate_jet_;
};

NLLS_REGISTER_REGISTER(VertexBase)
#define NLLS_REGISTER_VERTEX(name) \
    NLLS_REGISTER_CLASS(VertexBase, name)

template <class TupleType, int i>
struct VertexTypeAt {
    using VertexType = typename VertexTypeAt<typename TupleType::BaseType, i - 1>::VertexType;
};
template <class TupleType>
struct VertexTypeAt<TupleType, 0> {
    using VertexType = typename TupleType::VertexType;
};

template <class TupleType, int i>
struct VertexAt {
    static NLLS_INLINE typename VertexTypeAt<TupleType, i>::VertexType *Vertex(TupleType *t) {
        return VertexAt<typename TupleType::BaseType, i - 1>::Vertex(t);
    }
};
template <class TupleType>
struct VertexAt<TupleType, 0> {
    static NLLS_INLINE typename VertexTypeAt<TupleType, 0>::VertexType *Vertex(TupleType *t) {
        return t->vptr_;
    }
};
template <class ...V> class VertexTuple;
template <> class VertexTuple<> { template <class TupleType, int i> friend struct VertexAt; };
template <class V>
class VertexTuple<V> {
public:
    static constexpr int N = 1;
    using ThisType = VertexTuple<V>;
    using VertexType = V;
    explicit VertexTuple() { vptr_ = nullptr; }
    explicit VertexTuple(V *ptr) { vptr_ = ptr; }
    void SetVertices(V *ptr) { vptr_ = ptr; }
    template <int i> typename VertexTypeAt<ThisType, i>::VertexType *GetVertex() { return VertexAt<ThisType, i>::Vertex(this); }
    const VertexBase *GetVertex(int i) const { if (i == 0) { return vptr_; } else { return nullptr; } }
    VertexBase *GetVertex(int i) { if (i == 0) { return vptr_; } else { return nullptr; } }
    int GetVertexId(int i) const { if (i == 0) { return vptr_->GetId(); } else { -1; } }
    template <int i> static constexpr int PerturbationDim() { return (i == 0) ? V::DOF : 0; }
    template <int i> static constexpr int PerturbationOffset() { return (i == 0) ? 0 : V::DOF; }
    static constexpr int PerturbationDim(int i) { return (i == 0) ? V::DOF : 0; }
    static constexpr int PerturbationOffset(int i) { return (i == 0) ? 0 : V::DOF; }
    static constexpr int TotalPerturbationDim = PerturbationOffset<N>();
private:
    V *vptr_;
    template <class TupleType, int i> friend struct VertexAt;
};
template <typename V, typename ...VRest>
class VertexTuple<V, VRest...> : private VertexTuple<VRest...> {
public:
    static constexpr int N = sizeof...(VRest) + 1;
    using ThisType = VertexTuple<V, VRest...>;
    using BaseType = VertexTuple<VRest...>;
    using VertexType = V;
    VertexTuple() { vptr_ = nullptr; }
    VertexTuple(V *vptr, VRest*... vrest_ptrs) : BaseType(vrest_ptrs...) { vptr_ = vptr; }
    void SetVertices(V *vptr, VRest*... vrest_ptrs) { BaseType::SetVertices(vrest_ptrs...); vptr_ = vptr; }
    template <int i> typename VertexTypeAt<ThisType, i>::VertexType *GetVertex() { return VertexAt<ThisType, i>::Vertex(this); }
    const VertexBase *GetVertex(int i) const { if (i == 0) { return vptr_; } else { return BaseType::GetVertex(i - 1); } }
    VertexBase *GetVertex(int i) { if (i == 0) { return vptr_; } else { return BaseType::GetVertex(i - 1); } }
    int GetVertexId(int i) const { if (i == 0) { return vptr_->GetId(); } else { return BaseType::GetVertexId(i - 1); } }
    template <int i> static constexpr int PerturbationDim() { return (i == 0) ? V::DOF : BaseType::template PerturbationDim<i - 1>(); }
    template <int i> static constexpr int PerturbationOffset() { return (i == 0) ? 0 : BaseType::template PerturbationOffset<i - 1>() + V::DOF; }
    static constexpr int PerturbationDim(int i) { return (i == 0) ? V::DOF : BaseType::PerturbationDim(i - 1); }
    static constexpr int PerturbationOffset(int i) { return (i == 0) ? 0 : BaseType::PerturbationOffset(i - 1) + V::DOF; }
    static constexpr int TotalPerturbationDim = PerturbationOffset<N>();
private:
    V *vptr_;
    template <typename TupleType, int i> friend struct VertexAt;
};

} // namespace nlls