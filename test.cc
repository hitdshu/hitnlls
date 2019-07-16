#include "matrix/matrixx.h"
#include "matrix/matrix.h"
#include "matrix/vecxs.h"
#include "matrix/matrixs.h"
#include <chrono>
#include "node/node_point2f.h"
#include "node/node_rot2.h"
#include "node/node_pose2.h"
#include "factor/factor_se2_se2.h"
#include "factor/factor_point2f_se2.h"
#include "graph/naive_graph.h"
#include "solver/gn_solver.h"

using namespace ::hitnlls::matrix;
using namespace ::hitnlls::node;
using namespace ::hitnlls::factor;
using namespace ::hitnlls::graph;
using namespace ::hitnlls::solver;

namespace {
class Timer {
public:
    void Tic() {
        start_ticking_ = true;
        start_ = std::chrono::high_resolution_clock::now();
    }

    double Toc(std::string desc="") {
        if(!start_ticking_)
            return 0;
        start_ticking_ = false;
        end_ = std::chrono::high_resolution_clock::now();
        double t = std::chrono::duration<double, std::milli>(end_ - start_).count();
        std::cout << "Timer: " << t << " ms in " << desc << std::endl;
        return t;
    }

private:
    bool start_ticking_ = false;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_;
};
}

void TestMatrixxf() {
    Matrix33f m1;
    m1(0, 0) = 8.8902;
    m1(0, 1) = -2.6261;
    m1(0, 2) = 4.5115;
    m1(1, 0) = -2.6261;
    m1(1, 1) = 5.4043;
    m1(1, 2) = -7.4301;
    m1(2, 0) = 4.5115;
    m1(2, 1) = -7.4301;
    m1(2, 2) = 28.8281;
    Vector3f b;
    b[0] = 1;
    b[1] = 0;
    b[2] = 0;
    Vector3f x;
    x = m1 * b;
    (m1.SolveAxb(b)).Print("result");
    (m1.Inverse()).Print("inverse");
    (m1.Inverse().Transpose()).Print("inverse transpose");
    Matrix33f m2 = Matrix33f::Identity();
    m2.Print("Identity");
    (5 * m2).Print();
}

void TestVecxs() {
    Vecxs<Matrixxf> vxf1(100);
    Vecxs<Matrixxf> vxf2(100);
    Matrix33f tmp1;
    Matrix33f tmp2;
    vxf1[0] = tmp1 + 5;
    vxf2[99] = tmp2 + 1000;
    vxf2[0] = Matrix33f::Identity() * 6;
    
    Vecxs<Matrixxf> vxf(100);
    vxf = vxf1 - vxf2;
    vxf.Print();
}

void TestMatrixs() {
    Matrixsf ms1(3, 3);
    ms1(0, 0) = 0.05;
    ms1(0, 1) = 0.02;
    ms1(1, 0) = 0.02;
    ms1(1, 1) = 0.01;
    ms1(2, 2) = 1;
    Matrixsf ms2(3, 4);
    ms2(0, 0) = 2;
    ms2(1, 2) = 5;
    ms2(2, 2) = 9;
    Matrixsf ms = ms1 * ms2;
    ms1.SetLambdalm(2);
    Matrixsf l = ms1.CholeskyLLT();
    ms1.Print("original");
    l.Print("cholesky");
}

void TestMatrixsChol() {
    Matrixsxf sf(2, 2);

    Matrixxf mf1(2, 2);
    mf1(0, 0) = 0.1;
    mf1(0, 1) = 0.02;
    mf1(1, 0) = 0.02;
    mf1(1, 1) = 0.3;
    sf(0, 0) = mf1;
    Matrixxf mf2(2, 2);
    mf2(0, 0) = 0.4;
    mf2(0, 1) = 0.2;
    mf2(1, 0) = 0.2;
    mf2(1, 1) = 0.4;
    sf(1, 1) = mf2;

    Matrixsxf chol = sf.CholeskyLLT();

    sf.Print("original");
    chol.Print("cholesky");
}

void TestMatrixsCholHard() {
    Matrixsxf sf(4, 4);

    Matrixxf b00(2, 2);
    b00(0, 0) = 33.4328;
    b00(0, 1) = 11.1395;
    b00(1, 0) = 11.1395;
    b00(1, 1) = 30.0546;

    Matrixxf b10(3, 2);
    b10(0, 0) = -0.8622;
    b10(0, 1) = -1.3941;
    b10(1, 0) = 7.3386;
    b10(1, 1) = 13.4468;
    b10(2, 0) = -0.3757;
    b10(2, 1) = -1.6571;

    Matrixxf b11(3, 3);
    b11(0, 0) = 26.6782;
    b11(0, 1) = -6.3518;
    b11(0, 2) = 2.9747;
    b11(1, 0) = -6.3518;
    b11(1, 1) = 26.1955;
    b11(1, 2) = 5.5299;
    b11(2, 0) = 2.9747;
    b11(2, 1) = 5.5299;
    b11(2, 2) = 17.4060;

    Matrixxf b22(2, 2);
    b22(0, 0) = 23.6090;
    b22(0, 1) = -3.7071;
    b22(1, 0) = -3.7071;
    b22(1, 1) = 14.7363;

    Matrixxf b30(1, 2);
    b30(0, 0) = -0.0810;
    b30(0, 1) = 6.4773;

    Matrixxf b33(1, 1);
    b33(0, 0) = 14.3292;

    for (int idx = 0; idx < 1; ++idx) {
        sf(idx * 4 + 0, idx * 4 + 0) = b00;
        sf(idx * 4 + 1, idx * 4 + 0) = b10;
        sf(idx * 4 + 0, idx * 4 + 1) = b10.Transpose();
        sf(idx * 4 + 1, idx * 4 + 1) = b11;
        sf(idx * 4 + 2, idx * 4 + 2) = b22;
        sf(idx * 4 + 3, idx * 4 + 0) = b30;
        sf(idx * 4 + 0, idx * 4 + 3) = b30.Transpose();
        sf(idx * 4 + 3, idx * 4 + 3) = b33;
    }

    // sf.Print("original");
    sf.SetLambdalm(2.0);

    Timer timer;
    timer.Tic();
    Matrixsxf chol;
    for (int idx = 0; idx < 1000; ++idx) {
        chol = sf.CholeskyLLT();
    }
    timer.Toc("1000 cholesky");
    chol.Print("cholesky");
}

void TestMatrixsSolve() {
    Matrixsf ms(3, 3);
    ms(0, 0) = 0.05;
    ms(0, 1) = 0.02;
    ms(1, 0) = 0.02;
    ms(1, 1) = 0.05;
    ms(2, 2) = 1;
    ms.SetLambdalm(2);
    Matrixsf l = ms.CholeskyLLT();
    ms.Print("original");
    l.Print("cholesky");
    Vecfs b(3);
    b[0] = 1.1;
    b[1] = 2.2;
    b[2] = 3.3;
    Vecfs x = ms.SolveWithlm(b);
    std::cout << "Solver result " << x << std::endl;
}

void TestMatrixsSolveHard() {
    Matrixsxf sf(4, 4);
    Matrixxf b00(2, 2);
    b00(0, 0) = 33.4328;
    b00(0, 1) = 11.1395;
    b00(1, 0) = 11.1395;
    b00(1, 1) = 30.0546;
    Matrixxf b10(3, 2);
    b10(0, 0) = -0.8622;
    b10(0, 1) = -1.3941;
    b10(1, 0) = 7.3386;
    b10(1, 1) = 13.4468;
    b10(2, 0) = -0.3757;
    b10(2, 1) = -1.6571;
    Matrixxf b11(3, 3);
    b11(0, 0) = 26.6782;
    b11(0, 1) = -6.3518;
    b11(0, 2) = 2.9747;
    b11(1, 0) = -6.3518;
    b11(1, 1) = 26.1955;
    b11(1, 2) = 5.5299;
    b11(2, 0) = 2.9747;
    b11(2, 1) = 5.5299;
    b11(2, 2) = 17.4060;
    Matrixxf b22(2, 2);
    b22(0, 0) = 23.6090;
    b22(0, 1) = -3.7071;
    b22(1, 0) = -3.7071;
    b22(1, 1) = 14.7363;
    Matrixxf b30(1, 2);
    b30(0, 0) = -0.0810;
    b30(0, 1) = 6.4773;
    Matrixxf b33(1, 1);
    b33(0, 0) = 14.3292;
    for (int idx = 0; idx < 1; ++idx) {
        sf(idx * 4 + 0, idx * 4 + 0) = b00;
        sf(idx * 4 + 1, idx * 4 + 0) = b10;
        sf(idx * 4 + 0, idx * 4 + 1) = b10.Transpose();
        sf(idx * 4 + 1, idx * 4 + 1) = b11;
        sf(idx * 4 + 2, idx * 4 + 2) = b22;
        sf(idx * 4 + 3, idx * 4 + 0) = b30;
        sf(idx * 4 + 0, idx * 4 + 3) = b30.Transpose();
        sf(idx * 4 + 3, idx * 4 + 3) = b33;
    }
    sf.SetLambdalm(2.0);

    Vecxsxf vf(4);
    Matrixxf v0(2);
    v0[0] = 1.1;
    v0[1] = 2.2;
    Matrixxf v1(3);
    v1[0] = 3.3;
    v1[1] = 4.4;
    v1[2] = 5.5;
    Matrixxf v2(2);
    v2[0] = 6.6;
    v2[1] = 7.7;
    Matrixxf v3(1);
    v3[0] = 8.8;
    for (int idx = 0; idx < 1; ++idx) {
        vf[idx * 4 + 0] = v0;
        vf[idx * 4 + 1] = v1;
        vf[idx * 4 + 2] = v2;
        vf[idx * 4 + 3] = v3;
    }

    Timer timer;
    Vecxsxf xr;
    timer.Tic();
    for (int idx = 0; idx < 1000; ++idx) {
        xr = sf.SolveWithlm(vf);
        // Matrixsxf chol = sf.CholeskyLLT();
    }
    timer.Toc("1000 solver");

    std::cout << "Solved x " << xr;
    // Matrix<float, 2, 1> vf21 = v0;
    // std::cout << "vf21 " << vf21 << std::endl;
}

void TestNodePoint2f() {
    NodePose2 pos2;
    std::cout << "Before update " << pos2.GetEstimate() << std::endl;
    Vector3f update;
    update[0] = 0.7854;
    update[1] = 2;
    update[2] = 1;
    pos2.UpdatePlus(update);
    std::cout << "After update " << pos2.GetEstimate() << std::endl;
}

void TestFactorDerivSE2SE2() {
    SE2 twc1(SO2(0.5236));
    SE2 tcw1 = twc1.Inverse();
    Vector3f tcw1_noise;
    // tcw1_noise[0] = 1e-3;
    SE2 tcw1_meas = tcw1 + tcw1_noise;
    
    Vector2f pc1c2;
    pc1c2[0] = 1.1;
    pc1c2[1] = 2.2;
    SE2 tc1c2(SO2(0.7854), pc1c2);
    SE2 twc2 = twc1 * tc1c2;
    SE2 tcw2 = twc2.Inverse();
    Vector3f tcw2_noise;
    tcw2_noise[0] = 1e-3;
    SE2 tcw2_meas = tcw2 + tcw2_noise;

    NodePose2 pcw1;
    pcw1.SetEstimate(tcw1_meas);
    NodePose2 pcw2;
    pcw2.SetEstimate(tcw2_meas);

    FactorSE2SE2 fc1c2;
    fc1c2.SetNode(0, &pcw1);
    fc1c2.SetNode(1, &pcw2);

    std::cout << "Evaluate: " << fc1c2.Evaluate() << std::endl;
    std::cout << "Jacobian tc1w: " << fc1c2.Jacobian(0) << std::endl;
    std::cout << "Jacobian tc2w: " << fc1c2.Jacobian(1) << std::endl;
}

void TestFactorDerivPoint2fSE2() {
    SE2 twc1(SO2(0.5236));
    SE2 tcw1 = twc1.Inverse();
    Vector2f pc1c2;
    pc1c2[0] = 1.1;
    pc1c2[1] = 2.2;
    SE2 tc1c2(SO2(0.7854), pc1c2);
    SE2 twc = twc1 * tc1c2;
    SE2 tcw = twc.Inverse();
    Vector3f tcw_noise;
    tcw_noise[2] = 1e-3;
    tcw = tcw + tcw_noise;
    NodePose2 pcw;
    pcw.SetEstimate(tcw);

    ::hitnlls::matrix::Vector2f vw;
    vw[0] = 8.1;
    vw[1] = 7.5;
    NodePoint2f pw;
    pw.SetEstimate(vw);

    FactorPoint2fSE2 fps;
    fps.SetNode(0, &pw);
    fps.SetNode(1, &pcw);

    std::cout << "Evaluate: " << fps.Evaluate() << std::endl;
    std::cout << "Jacobian pw: " << fps.Jacobian(0) << std::endl;
    std::cout << "Jacobian tcw: " << fps.Jacobian(1) << std::endl;
    std::cout << "Jacobian tcw subset: " << fps.Jacobian(1).Block(0, 0, 2, 1) << std::endl;
    std::cout << "Jacobian tcw subset: " << fps.Jacobian(1).SetBlock(0, 0, 2, 1, vw) << std::endl;
}

void TestNaiveGraph() {
    NaiveGraph ng;
}

int main(int argc, char **argv) {
    // TestMatrixxf();
    // TestVecxs();
    // TestMatrixs();
    // TestMatrixsChol();
    // TestMatrixsCholHard();
    // TestMatrixsSolve();
    // TestMatrixsSolveHard();
    // TestNodePoint2f();
    // TestFactorDerivSE2SE2();
    // TestFactorDerivPoint2fSE2();
    TestNaiveGraph();
}