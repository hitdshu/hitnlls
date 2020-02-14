#pragma once

#include <vector>
#include "factor/loss.h"
#include "factor/equality.h"
#include "factor/inequality.h"

namespace hitcadmm {
    
class Problem {
public:
    typedef std::shared_ptr<Problem> Ptr;

    void AddLoss(const Loss::Ptr &loss) { losses_.push_back(loss); }
    void AddEquality(const Equality::Ptr &equa) { equas_.push_back(equa); }
    void AddInequality(const Inequality::Ptr &ine) { ines_.push_back(ine); }

    void Build();
    void SetEstimates(const Eigen::VectorXd &est);

    Eigen::MatrixXd GetP() const { return P_; }
    Eigen::VectorXd Getq() const { return q_; }
    Eigen::MatrixXd GetA() const { return A_; }
    Eigen::VectorXd Getb() const { return b_; }
    Eigen::MatrixXd GetC() const { return C_; }
    Eigen::VectorXd Getd() const { return d_; }
    int GetXDim() const { return xdim_; }
    int GetEDim() const { return edim_; }
    int GetSDim() const { return sdim_; }
    std::vector<Vertex::Ptr> GetCones() const { return cones_; }

private:
    std::vector<Loss::Ptr> losses_;
    std::vector<Equality::Ptr> equas_;
    std::vector<Inequality::Ptr> ines_;
    std::vector<Vertex::Ptr> xs_;

    Eigen::MatrixXd P_;
    Eigen::VectorXd q_;
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    Eigen::MatrixXd C_;
    Eigen::VectorXd d_;
    int xdim_;
    int edim_;
    int sdim_;
    std::vector<Vertex::Ptr> cones_;
};

void Problem::Build() {
    std::vector<int> xi;
    std::vector<int> ai;
    std::vector<int> ci;
    xdim_ = 0;
    for (auto &loss : losses_) {
        Vertex::Ptr v = loss->GetV();
        xs_.push_back(v);
        xi.push_back(xdim_);
        xdim_ += v->GetDim();
    }
    edim_ = 0;
    for (auto &eq : equas_) {
        ai.push_back(edim_);
        edim_ += eq->Getb().size();
    }
    sdim_ = 0;
    for (auto &ie : ines_) {
        ci.push_back(sdim_);
        sdim_ += ie->Getd().size();
    }
    P_ = Eigen::MatrixXd(xdim_, xdim_);
    P_.setZero();
    q_ = Eigen::VectorXd(xdim_);
    q_.setZero();
    A_ = Eigen::MatrixXd(edim_, xdim_);
    A_.setZero();
    b_ = Eigen::VectorXd(edim_);
    b_.setZero();
    C_ = Eigen::MatrixXd(sdim_, xdim_);
    C_.setZero();
    d_ = Eigen::VectorXd(sdim_);
    d_.setZero();
    for (auto &loss : losses_) {
        Vertex::Ptr v = loss->GetV();
        int vidx = (std::find(xs_.begin(), xs_.end(), v) - xs_.begin());
        int tidx = xi[vidx];
        int tdim = v->GetDim();
        P_.block(tidx, tidx, tdim, tdim) = loss->GetP();
        q_.block(tidx, 0, tdim, 1) = loss->Getq();
    }
    int ax = 0;
    for (auto &eq : equas_) {
        Vertex::Ptr v = eq->GetV();
        int vidx = (std::find(xs_.begin(), xs_.end(), v) - xs_.begin());
        int xidx = xi[vidx];
        int aidx = ai[ax];
        int tdim = v->GetDim();
        A_.block(aidx, xidx, eq->Getb().size(), tdim) = eq->GetA();
        b_.block(aidx, 0, eq->Getb().size(), 1) = eq->Getb();
        ax++;
    }
    int cx = 0;
    cones_.clear();
    for (auto &ie : ines_) {
        Vertex::Ptr v = ie->GetV();
        int vidx = (std::find(xs_.begin(), xs_.end(), v) - xs_.begin());
        int xidx = xi[vidx];
        int cidx = ci[cx];
        int tdim = v->GetDim();
        C_.block(cidx, xidx, ie->Getd().size(), tdim) = ie->GetC();
        d_.block(cidx, 0, ie->Getd().size(), 1) = ie->Getd();
        cones_.push_back(ie->GetCone());
        cx++;
    }
}

void Problem::SetEstimates(const Eigen::VectorXd &est) {
    int idx = 0;
    for (auto &x : xs_) {
        int xdim = x->GetDim();
        x->SetVector(est.block(idx, 0, xdim, 1));
        idx += xdim;
    }
}

} // namespace hitcadmm