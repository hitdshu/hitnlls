#pragma once

#include "common/register.h"
#include "type/vertex.h"
#include "type/kernel.h"

namespace hitnlls {

class FactorBase {
public:
    FactorBase() { kernel_ = nullptr; }
    virtual ~FactorBase() = default;

    void SetKernel(KernelBase *kernel) { kernel_ = kernel; }
    KernelBase *GetKernel() const { return kernel_; }

    virtual int NumVertices() const = 0;
    virtual VertexBase *GetVertexAt(int i) = 0;
    virtual const VertexBase *GetVertexAt(int i) const = 0;
    virtual int GetVertexIdAt(int i) const = 0;

    virtual void Compute() = 0;
    virtual double GetChi() = 0;
    virtual matrix::MatrixXd JacobTInfoRes(int i) = 0;
    virtual matrix::MatrixXd JacobTInfoJacob(int i, int j) = 0;

    FactorBase(const FactorBase &) = delete;
    FactorBase &operator=(const FactorBase &) = delete;
    FactorBase &operator=(const FactorBase &&) = delete;

protected:
    void Robustify();

    matrix::Vector3d robust_scales_;

    KernelBase *kernel_;
};

template <int EDIM, typename ...V>
class Factor : public FactorBase {
public:
    using VertexTupleType = VertexTuple<V...>;
    using ThisType = Factor<EDIM, V...>;
    using BaseType = FactorBase;
    static constexpr int VertexNums = VertexTupleType::VertexNums;
    static constexpr int TotalPerturbationDim = VertexTupleType::TotalPerturbationDim;
    static constexpr int ErrorDim = EDIM;
    using ErrorVectorType = matrix::Matrix<double, EDIM, 1>;
    using InformationMatrixType = matrix::Matrix<double, EDIM, EDIM>;
    using JacobianMatrixType = matrix::Matrix<double, EDIM, TotalPerturbationDim>;

    Factor() { information_.SetIdentity(); }
    Factor(V*... vptrs) : vertices_(vptrs...) { information_.SetIdentity(); }

    void SetVertices(V*... vptrs) { vertices_.SetVertices(vptrs...); }
    void SetInformation(const InformationMatrixType &info) { information_ = info; }
    void SetMeasurement(const ErrorVectorType &meas) { measure_ = meas; }
    const VertexTupleType &GetVertexTuple() const { return vertices_; }
    VertexTupleType &GetVertexTuple() { return vertices_; }

    virtual void ErrorAndJacobian(bool chi_only = false) = 0;

    virtual int NumVertices() const override final { return VertexNums; }
    virtual VertexBase *GetVertexAt(int i) override final { return vertices_.GetVertex(i); }
    virtual const VertexBase *GetVertexAt(int i) const override final { return vertices_.GetVertex(i); }
    virtual int GetVertexIdAt(int i) const override final { return vertices_.GetVertexId(i); }

    virtual void Compute() { ErrorAndJacobian(); BaseType::Robustify(); }
    virtual double GetChi() override final { return chi_; }
    virtual matrix::MatrixXd JacobTInfoRes(int i) override final {
        return jacobians_.Block(0, vertices_.PerturbationOffset(i), EDIM, vertices_.PerturbationOffset(i)).Transpose() * information_
            * residual_ * BaseType::robust_scales_[1];
    }
    virtual matrix::MatrixXd JacobTInfoJacob(int i, int j) override final {
        jacobians_.Block(0, vertices_.PerturbationOffset(i), EDIM, vertices_.PerturbationOffset(i)).Transpose() * information_
            * jacobians_.Block(0, vertices_.PerturbationOffset(j), EDIM, vertices_.PerturbationOffset(j)) * BaseType::robust_scales_[1];
    }

public:
    VertexTupleType vertices_;
    InformationMatrixType information_;
    ErrorVectorType measure_;
    ErrorVectorType residual_;
    JacobianMatrixType jacobians_;
    double chi_;
};

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
    static constexpr int VertexNums = VertexTupleType::VertexNums;
    static constexpr int TotalPerturbationDim = VertexTupleType::TotalPerturbationDim;
    static constexpr int ErrorDim = EDIM;
    using ErrorVectorTypeJet = matrix::Matrix<geometry::Jetd, EDIM, 1>;
    using EvalVectorTypeJet = ErrorVectorTypeJet;

    virtual EvalVectorTypeJet operator()(VertexTupleType &vars) = 0;

    template <class ThisType, int i>
    friend struct JacobianUpdater;

    virtual void ErrorAndJacobian(bool chi_only = false) override final {
        EvalVectorTypeJet eval_ad = this->operator()(BaseType::GetVertexTuple());
        geometry::ConvertMatrix(BaseType::residual_, eval_ad);
        BaseType::residual_ = BaseType::measure_ - BaseType::residual_;
        auto tmp = BaseType::residual_.Transpose() * BaseType::information_ * BaseType::residual_;
        BaseType::chi_ = tmp(0, 0);
        if (chi_only) {
            return;
        }
        JacobianUpdater<ThisType, VertexNums - 1>::Update(*this);
    }

    template <int i>
    void JacobianUpdateAt() {
        constexpr int PertDim = VertexTupleType::template PerturbationDim<i>();
        using VertexPerturbation = matrix::Matrix<geometry::Jetd, PertDim, 1>;
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
    static constexpr int VertexNums = VertexTupleType::VertexNums;
    static constexpr int TotalPerturbationDim = VertexTupleType::TotalPerturbationDim;
    static constexpr int ErrorDim = EDIM;
    using ErrorVectorType = typename BaseType::ErrorVectorType;
    using EvalVectorType = ErrorVectorType;

    virtual EvalVectorType operator()(VertexTupleType &vars) = 0;

    template <class ThisType, int i>
    friend struct JacobianUpdater;

    virtual void ErrorAndJacobian(bool chi_only = false) override final {
        EvalVectorType eval = this->operator()(BaseType::GetVertexTuple());
        BaseType::residual_ = BaseType::measure_ - eval;
        auto tmp = BaseType::residual_.Transpose() * BaseType::information_ * BaseType::residual_;
        BaseType::chi_ = tmp(0, 0);
        if (chi_only) {
            return;
        }
        JacobianUpdater<ThisType, VertexNums - 1>::Update(*this);
    }

    template <int i>
    void JacobianUpdateAt() {
        constexpr int PertDim = VertexTupleType::template PerturbationDim<i>();
        using VertexPerturbation = matrix::Matrix<double, PertDim, 1>;
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

HITNLLS_REGISTER_REGISTER(FactorBase)
#define HITNLLS_REGISTER_FACTOR(name) \
    HITNLLS_REGISTER_CLASS(FactorBase, name)

} // namespace hitnlls