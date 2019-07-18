#include "catch/catch.hpp"
#include "matrix/matrix.h"

TEST_CASE( "Test dense matrix operation", "[matrix]") {
    ::hitnlls::matrix::Matrix33f mat;
    mat(0, 0) = 8.8902;
    mat(0, 1) = -2.6261;
    mat(0, 2) = 4.5115;
    mat(1, 0) = -2.6261;
    mat(1, 1) = 5.4043;
    mat(1, 2) = -7.4301;
    mat(2, 0) = 4.5115;
    mat(2, 1) = -7.4301;
    mat(2, 2) = 28.8281;

    ::hitnlls::matrix::Vector3f b;
    b[0] = 1;
    b[1] = 0;
    b[2] = 0;

    ::hitnlls::matrix::Vector3f x = mat.SolveAxb(b);
    REQUIRE(x[0] == Approx(0.132092));
    REQUIRE(x[1] == Approx(0.0553959));
    REQUIRE(x[2] == Approx(-0.00639428));

    ::hitnlls::matrix::Matrix33f mat_inv = mat.Inverse();
    REQUIRE(mat_inv(0, 0) == Approx(0.132092));
    REQUIRE(mat_inv(0, 1) == Approx(0.0553959));
    REQUIRE(mat_inv(0, 2) == Approx(-0.00639428));
    REQUIRE(mat_inv(1, 0) == Approx(0.0553959));
    REQUIRE(mat_inv(1, 1) == Approx(0.309824));
    REQUIRE(mat_inv(1, 2) == Approx(0.0711841));
    REQUIRE(mat_inv(2, 0) == Approx(-0.00639428));
    REQUIRE(mat_inv(2, 1) == Approx(0.0711841));
    REQUIRE(mat_inv(2, 2) == Approx(0.0540359));

    ::hitnlls::matrix::Matrix33f mat_inv_tran = mat_inv.Transpose();
    REQUIRE(mat_inv_tran(0, 0) == Approx(0.132092));
    REQUIRE(mat_inv_tran(0, 1) == Approx(0.0553959));
    REQUIRE(mat_inv_tran(0, 2) == Approx(-0.00639428));
    REQUIRE(mat_inv_tran(1, 0) == Approx(0.0553959));
    REQUIRE(mat_inv_tran(1, 1) == Approx(0.309824));
    REQUIRE(mat_inv_tran(1, 2) == Approx(0.0711841));
    REQUIRE(mat_inv_tran(2, 0) == Approx(-0.00639428));
    REQUIRE(mat_inv_tran(2, 1) == Approx(0.0711841));
    REQUIRE(mat_inv_tran(2, 2) == Approx(0.0540359));

    ::hitnlls::matrix::Matrix33f mat_diag5 = 5 * ::hitnlls::matrix::Matrix33f::Identity();
    REQUIRE(mat_diag5(0, 0) == Approx(5));
    REQUIRE(mat_diag5(0, 1) == Approx(0));
    REQUIRE(mat_diag5(0, 2) == Approx(0));
    REQUIRE(mat_diag5(1, 0) == Approx(0));
    REQUIRE(mat_diag5(1, 1) == Approx(5));
    REQUIRE(mat_diag5(1, 2) == Approx(0));
    REQUIRE(mat_diag5(2, 0) == Approx(0));
    REQUIRE(mat_diag5(2, 1) == Approx(0));
    REQUIRE(mat_diag5(2, 2) == Approx(5));
}