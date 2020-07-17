#pragma once

#include "../common/register.h"

#include "vertex.h"
#include "kernel.h"

namespace nlls {

class FactorBase {
public:
    NLLS_NONCOPYABLE(FactorBase)
    enum Status { Active, NonActive };
    FactorBase() { kernel_ = nullptr; status_= Active; scale_ = 1; }
    virtual ~FactorBase() = default;

    void SetStatus(Status st) { status_ = st; }
    Status GetStatus() const { return status_; }
    void SetKernel(KernelBase *kernel) { kernel_ = kernel; }
    KernelBase *GetKernel() const { return kernel_; }
    float GetChi() const { return chi_; }
    MatrixXf JacobTransRes(int i) const { return Jacobian(i).Transpose() * Residual(); }
    MatrixXf JacobTransJacob(int i, int j) const { MatrixXf ji = Jacobian(i); MatrixXf jj = Jacobian(j); return ji.Transpose() * jj; }

    virtual int NumObservations() const = 0;
    virtual int NumVertices() const = 0;
    virtual VertexBase *GetVertexAt(int i) = 0;
    virtual const VertexBase *GetVertexAt(int i) const = 0;
    virtual int GetVertexIdAt(int i) const = 0;
    virtual int GetVertexDimOffset(int i) const = 0;
    virtual void Compute(bool chi_only = false) = 0;
    virtual MatrixXf Residual() const = 0;
    virtual MatrixXf Jacobian(int i) const = 0;
    virtual MatrixXf Jacobian() const = 0;

protected:
    void Robustify();

    float chi_;
    float scale_;
    Status status_;
    KernelBase *kernel_;
};

template <int EDIM, typename ...V>
class Factor : public FactorBase {
public:
    using VertexTupleType = VertexTuple<V...>;
    using ThisType = Factor<EDIM, V...>;
    using BaseType = FactorBase;
    static constexpr int VertexNums = VertexTupleType::N;
    static constexpr int TotalPerturbationDim = VertexTupleType::TotalPerturbationDim;
    static constexpr int ErrorDim = EDIM;
    using ErrorVector = Matrix<float, EDIM, 1>;
    using SqrtInfoMatrix = Matrix<float, EDIM, EDIM>;
    using JacobianMatrix = Matrix<float, EDIM, TotalPerturbationDim>;

    Factor() { sqrt_info_.SetIdentity(); }
    Factor(V*... vptrs) : vertices_(vptrs...) { sqrt_info_.SetIdentity(); }

    void SetVertices(V*... vptrs) { vertices_.SetVertices(vptrs...); }
    void SetSqrtInfo(const SqrtInfoMatrix &sqrt_info) { sqrt_info_ = sqrt_info; }
    void SetMeasurement(const ErrorVector &meas) { measure_ = meas; }
    const VertexTupleType &GetVertexTuple() const { return vertices_; }
    VertexTupleType &GetVertexTuple() { return vertices_; }

    virtual int NumObservations() const override final { return EDIM; }
    virtual int NumVertices() const override final { return VertexNums; }
    virtual VertexBase *GetVertexAt(int i) override final { return vertices_.GetVertex(i); }
    virtual const VertexBase *GetVertexAt(int i) const override final { return vertices_.GetVertex(i); }
    virtual int GetVertexIdAt(int i) const override final { return vertices_.GetVertexId(i); }
    virtual int GetVertexDimOffset(int i) const override final { return vertices_.PerturbationOffset(i); }
    virtual void Compute(bool chi_only = false) { ErrorAndJacobian(chi_only); BaseType::Robustify(); }
    virtual MatrixXf Residual() const override final { return residual_ * BaseType::scale_; }
    virtual MatrixXf Jacobian(int i) const override final {
        return jacobians_.Block(0, vertices_.PerturbationOffset(i), EDIM, vertices_.PerturbationDim(i)) * BaseType::scale_;
    }
    virtual MatrixXf Jacobian() const override final { return jacobians_ * BaseType::scale_; }

    virtual void ErrorAndJacobian(bool chi_only = false) = 0;

protected:
    VertexTupleType vertices_;
    SqrtInfoMatrix sqrt_info_;
    ErrorVector measure_;
    ErrorVector residual_;
    JacobianMatrix jacobians_;
};

NLLS_REGISTER_REGISTER(FactorBase)
#define NLLS_REGISTER_FACTOR(name) \
    NLLS_REGISTER_CLASS(FactorBase, name)

template <typename ErrorFactorType, int i>
struct JacobianUpdater {
    static inline void Update(ErrorFactorType& ef) {
        JacobianUpdater<ErrorFactorType, i - 1>::Update(ef);
        ef.template JacobianUpdateAt<i>();
    }
};
template <typename ErrorFactorType>
struct JacobianUpdater<ErrorFactorType, 0> {
    static inline void Update(ErrorFactorType& ef) {
        ef.template JacobianUpdateAt<0>();
    }
};

template <int EDIM, typename ...V>
class FactorAutoDiff : public Factor<EDIM, V...> {
public:
    using ThisType = FactorAutoDiff<EDIM, V...>;
    using BaseType = Factor<EDIM, V...>;
    using VertexTupleType = typename BaseType::VertexTupleType;
    static constexpr int VertexNums = VertexTupleType::N;
    static constexpr int TotalPerturbationDim = VertexTupleType::TotalPerturbationDim;
    static constexpr int ErrorDim = EDIM;
    using ErrorVectorJet = Matrix<Jetf, EDIM, 1>;
    using EvalVectorTypeJet = ErrorVectorJet;

    virtual EvalVectorTypeJet operator()(VertexTupleType &vars) = 0;

    template <class ErrorFactorType, int i>
    friend struct JacobianUpdater;

    virtual void ErrorAndJacobian(bool chi_only = false) override final {
        EvalVectorTypeJet eval_ad = this->operator()(BaseType::GetVertexTuple());
        ConvertMatrix(BaseType::residual_, eval_ad);
        BaseType::residual_ = BaseType::sqrt_info_ * (BaseType::measure_ - BaseType::residual_);
        BaseType::chi_ = pow(BaseType::residual_.Norm(), 2);
        if (chi_only) {
            return;
        }
        JacobianUpdater<ThisType, VertexNums - 1>::Update(*this);
        BaseType::jacobians_ = BaseType::sqrt_info_ * BaseType::jacobians_;
    }

    template <int i>
    void JacobianUpdateAt() {
        constexpr int PertDim = VertexTupleType::template PerturbationDim<i>();
        using VertexPerturbation = Matrix<Jetf, PertDim, 1>;
        VertexPerturbation pert;
        pert.SetZero();
        auto *var = BaseType::vertices_.template GetVertex<i>();
        var->SetEstimate(var->GetEstimate());
        const int jc_offset = BaseType::vertices_.template PerturbationOffset<i>();
        for (int c = 0; c < PertDim; ++c) {
            pert[c].Deriv() = 1.0;
            var->ApplyPerturbation(pert);
            EvalVectorTypeJet eval_ad = this->operator()(BaseType::GetVertexTuple());
            for (int r = 0; r < ErrorDim; ++r) {
                BaseType::jacobians_(r, c + jc_offset) = eval_ad[r].Deriv();
            }
            var->SetEstimate(var->GetEstimate());
            pert[c].Deriv() = 0;
        }
    }
};

template <int EDIM, typename ...V>
class FactorNumeDiff : Factor<EDIM, V...> {
public:
    using ThisType = FactorNumeDiff<EDIM, V...>;
    using BaseType = Factor<EDIM, V...>;
    using VertexTupleType = typename BaseType::VertexTupleType;
    static constexpr int VertexNums = VertexTupleType::N;
    static constexpr int TotalPerturbationDim = VertexTupleType::TotalPerturbationDim;
    static constexpr int ErrorDim = EDIM;
    using ErrorVector = typename BaseType::ErrorVector;
    using EvalVectorType = ErrorVector;

    virtual EvalVectorType operator()(VertexTupleType &vars) = 0;

    template <class ThisType, int i>
    friend struct JacobianUpdater;

    virtual void ErrorAndJacobian(bool chi_only = false) override final {
        EvalVectorType eval = this->operator()(BaseType::GetVertexTuple());
        BaseType::residual_ = BaseType::sqrt_info_ * (BaseType::measure_ - eval);
        BaseType::chi_ = pow(BaseType::residual_.Norm(), 2);
        if (chi_only) {
            return;
        }
        JacobianUpdater<ThisType, VertexNums - 1>::Update(*this);
        BaseType::jacobians_ = BaseType::sqrt_info_ * BaseType::jacobians_;
    }

    template <int i>
    void JacobianUpdateAt() {
        constexpr int PertDim = VertexTupleType::template PerturbationDim<i>();
        using VertexPerturbation = Matrix<float, PertDim, 1>;
        VertexPerturbation pert;
        pert.SetZero();
        auto *var = BaseType::vertices_.template GetVertex<i>();
        auto ori_est = var->GetEstimate();
        const int jc_offset = BaseType::vertices_.template PerturbationOffset<i>();
        for (int c = 0; c < PertDim; ++c) {
            pert[c] = -1e-7;
            var->ApplyPerturbation(pert);
            EvalVectorType eval1 = this->operator()(BaseType::GetVertexTuple());
            var->SetEstimate(ori_est);
            pert[c] = 1e-7;
            var->ApplyPerturbation(pert);
            EvalVectorType eval2 = this->operator()(BaseType::GetVertexTuple());
            var->SetEstimate(ori_est);
            for (int r = 0; r < ErrorDim; ++r) {
                BaseType::jacobians_(r, c + jc_offset) = (eval2[r] - eval1[r]) / (2e-7);
            }
            pert[c] = 0;
        }
    }
};

} // namespace nlls