#pragma once

#include <vector>

#include "admm/loss.h"
#include "admm/equality.h"
#include "admm/inequality.h"

namespace hitcadmm {
    
class Problem {
public:
    typedef std::shared_ptr<Problem> Ptr;

    void AddLoss(const Loss::Ptr &loss) { losses_.push_back(loss); }
    void AddEquality(const Equality::Ptr &equa) { equas_.push_back(equa); }
    void AddInequality(const Inequality::Ptr &ine) { ines_.push_back(ine); }

    void Build();
    void SetEstimates(const hitnlls::matrix::VectorXd &est);

    hitnlls::matrix::MatrixXd GetP() const { return P_; }
    hitnlls::matrix::VectorXd Getq() const { return q_; }
    hitnlls::matrix::MatrixXd GetA() const { return A_; }
    hitnlls::matrix::VectorXd Getb() const { return b_; }
    hitnlls::matrix::MatrixXd GetC() const { return C_; }
    hitnlls::matrix::VectorXd Getd() const { return d_; }
    int GetXDim() const { return xdim_; }
    int GetEDim() const { return edim_; }
    int GetSDim() const { return sdim_; }
    std::vector<Variable::Ptr> GetCones() const { return cones_; }

private:
    std::vector<Loss::Ptr> losses_;
    std::vector<Equality::Ptr> equas_;
    std::vector<Inequality::Ptr> ines_;
    std::vector<Variable::Ptr> xs_;

    hitnlls::matrix::MatrixXd P_;
    hitnlls::matrix::VectorXd q_;
    hitnlls::matrix::MatrixXd A_;
    hitnlls::matrix::VectorXd b_;
    hitnlls::matrix::MatrixXd C_;
    hitnlls::matrix::VectorXd d_;
    int xdim_;
    int edim_;
    int sdim_;
    std::vector<Variable::Ptr> cones_;
};

void Problem::Build() {
    std::vector<int> xi;
    std::vector<int> ai;
    std::vector<int> ci;
    xdim_ = 0;
    for (auto &loss : losses_) {
        Variable::Ptr v = loss->GetV();
        xs_.push_back(v);
        xi.push_back(xdim_);
        xdim_ += v->GetDim();
    }
    edim_ = 0;
    for (auto &eq : equas_) {
        ai.push_back(edim_);
        edim_ += eq->Getb().Size();
    }
    sdim_ = 0;
    for (auto &ie : ines_) {
        ci.push_back(sdim_);
        sdim_ += ie->Getd().Size();
    }
    P_ = hitnlls::matrix::MatrixXd(xdim_, xdim_);
    P_.SetZero();
    q_ = hitnlls::matrix::VectorXd(xdim_);
    q_.SetZero();
    A_ = hitnlls::matrix::MatrixXd(edim_, xdim_);
    A_.SetZero();
    b_ = hitnlls::matrix::VectorXd(edim_);
    b_.SetZero();
    C_ = hitnlls::matrix::MatrixXd(sdim_, xdim_);
    C_.SetZero();
    d_ = hitnlls::matrix::VectorXd(sdim_);
    d_.SetZero();
    for (auto &loss : losses_) {
        Variable::Ptr v = loss->GetV();
        int vidx = (std::find(xs_.begin(), xs_.end(), v) - xs_.begin());
        int tidx = xi[vidx];
        int tdim = v->GetDim();
        P_.Block(tidx, tidx, tdim, tdim) = loss->GetP();
        q_.Block(tidx, 0, tdim, 1) = loss->Getq();
    }
    int ax = 0;
    for (auto &eq : equas_) {
        Variable::Ptr v = eq->GetV();
        int vidx = (std::find(xs_.begin(), xs_.end(), v) - xs_.begin());
        int xidx = xi[vidx];
        int aidx = ai[ax];
        int tdim = v->GetDim();
        A_.Block(aidx, xidx, eq->Getb().Size(), tdim) = eq->GetA();
        b_.Block(aidx, 0, eq->Getb().Size(), 1) = eq->Getb();
        ax++;
    }
    int cx = 0;
    cones_.clear();
    for (auto &ie : ines_) {
        Variable::Ptr v = ie->GetV();
        int vidx = (std::find(xs_.begin(), xs_.end(), v) - xs_.begin());
        int xidx = xi[vidx];
        int cidx = ci[cx];
        int tdim = v->GetDim();
        C_.Block(cidx, xidx, ie->Getd().Size(), tdim) = ie->GetC();
        d_.Block(cidx, 0, ie->Getd().Size(), 1) = ie->Getd();
        cones_.push_back(ie->GetCone());
        cx++;
    }
}

void Problem::SetEstimates(const hitnlls::matrix::VectorXd &est) {
    int idx = 0;
    for (auto &x : xs_) {
        int xdim = x->GetDim();
        x->SetVector(est.Block(idx, 0, xdim, 1));
        idx += xdim;
    }
}

} // namespace hitcadmm