#include "catch/catch.hpp"
#include "matrix/matrixs.h"

TEST_CASE( "Test sparse matrix operation", "[matrixs]") {
    ::hitnlls::matrix::Matrixsf ms1(3, 3);
    ms1(0, 0) = 0.05;
    ms1(0, 1) = 0.02;
    ms1(1, 0) = 0.02;
    ms1(1, 1) = 0.01;
    ms1(2, 2) = 1;
    ::hitnlls::matrix::Matrixsf ms2(3, 4);
    ms2(0, 0) = 2;
    ms2(1, 2) = 5;
    ms2(2, 2) = 9;
    ::hitnlls::matrix::Matrixsf ms = ms1 * ms2;
    REQUIRE(ms(0, 0) == Approx(0.2));
    REQUIRE(ms(1, 0) == Approx(0.09));
    REQUIRE(ms(2, 0) == Approx(9));


    ms1.SetLambdalm(2);
    ::hitnlls::matrix::Matrixsf l = ms1.CholeskyLLT();
    REQUIRE(l(0, 0) == Approx(1.43178));
    REQUIRE(l(1, 0) == Approx(0.0139686));
    REQUIRE(l(1, 1) == Approx(1.41768));
    REQUIRE(l(2, 2) == Approx(1.73205));
}