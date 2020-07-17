#include "question.h"
#include "admm.h"

namespace nlls {

Question::~Question() {
    for (auto loss : losses_) {
        delete loss;
    }
    for (auto eq : equas_) {
        delete eq;
    }
    for (auto ie : ines_) {
        delete ie;
    }
    for (auto v : xs_) {
        delete v;
    }
    for (auto c : cones_) {
        delete c;
    }
}

void Question::Solve() {
    internal::Admm solver;
    solver.SetQuestion(this);
    solver.Solve();
}

void Question::Build() {
    std::vector<int> xi;
    std::vector<int> ai;
    std::vector<int> ci;
    xdim_ = 0;
    for (auto &loss : losses_) {
        Variable *v = loss->GetV();
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
    P_ = MatrixXf(xdim_, xdim_);
    P_.SetZero();
    q_ = VectorXf(xdim_);
    q_.SetZero();
    A_ = MatrixXf(edim_, xdim_);
    A_.SetZero();
    b_ = VectorXf(edim_);
    b_.SetZero();
    C_ = MatrixXf(sdim_, xdim_);
    C_.SetZero();
    d_ = VectorXf(sdim_);
    d_.SetZero();
    for (auto &loss : losses_) {
        using std::find;
        Variable *v = loss->GetV();
        int vidx = find(xs_.begin(), xs_.end(), v) - xs_.begin();
        int tidx = xi[vidx];
        int tdim = v->GetDim();
        P_.Block(tidx, tidx, tdim, tdim) = loss->GetP();
        q_.Block(tidx, 0, tdim, 1) = loss->Getq();
    }
    int ax = 0;
    for (auto &eq : equas_) {
        using std::find;
        Variable *v = eq->GetV();
        int vidx = find(xs_.begin(), xs_.end(), v) - xs_.begin();
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
        using std::find;
        Variable *v = ie->GetV();
        int vidx = find(xs_.begin(), xs_.end(), v) - xs_.begin();
        int xidx = xi[vidx];
        int cidx = ci[cx];
        int tdim = v->GetDim();
        C_.Block(cidx, xidx, ie->Getd().Size(), tdim) = ie->GetC();
        d_.Block(cidx, 0, ie->Getd().Size(), 1) = ie->Getd();
        cones_.push_back(ie->GetCone());
        cx++;
    }
}

void Question::SetEstimates(const VectorXf &est) {
    int idx = 0;
    for (auto &x : xs_) {
        int xdim = x->GetDim();
        x->SetVector(est.Block(idx, 0, xdim, 1));
        idx += xdim;
    }
}

} // namespace nlls