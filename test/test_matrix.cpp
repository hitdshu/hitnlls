#include "../matrix/dense.h"

#include "catch.hpp"

TEST_CASE("Matrix addition", "[matrix]") {
    nlls::Matrix33f m1 = nlls::Matrix33f::Random();
    nlls::Matrix33f m2 = nlls::Matrix33f::Random();
    nlls::Matrix33f m3 = m1 + m2;
    for (int r = 0; r < m3.Rows(); ++r) {
        for (int c = 0; c < m3.Cols(); ++c) {
            REQUIRE(m3(r, c) == Approx(m1(r, c) + m2(r, c)));
        }
    }
}

TEST_CASE("Matrix transpose and multiply", "[matrix]") {
    nlls::Matrix33d A;
    A << 0.537667139546100, 0.862173320368121, -0.433592022305684, 
        1.833885014595086, 0.318765239858981, 0.342624466538650, 
        -2.258846861003648, -1.307688296305273, 3.578396939725760;
    nlls::Matrix33d B;
    B << 2.769437029884877, 0.725404224946106, -0.204966058299775, 
        -1.349886940156521, -0.063054873189656, -0.124144348216312, 
        3.034923466331855, 0.714742903826096, 1.489697607785465;
    nlls::Matrix33d C = A.Transpose() * B;
    REQUIRE(C(0, 0) == Approx(-7.841929490249430));
    REQUIRE(C(0, 1) == Approx(-1.340104137130485));
    REQUIRE(C(0, 2) == Approx(-3.702868739301489));
    REQUIRE(C(1, 0) == Approx(-2.011296211759626));
    REQUIRE(C(1, 1) == Approx(-0.329336462746363));
    REQUIRE(C(1, 2) == Approx(-2.164349296718408));
    REQUIRE(C(2, 0) == Approx(9.196850749029156));
    REQUIRE(C(2, 1) == Approx(2.221500192569268));
    REQUIRE(C(2, 2) == Approx(5.377066117457115));
}

TEST_CASE("Matrix svd", "[matrix]") {
    nlls::Matrix43d A;
    A << 0.5377, 0.3188, 3.5784, 
        1.8339, -1.3077, 2.7694, 
        -2.2588, -0.4336, -1.3499, 
        0.8622, 0.3426, 3.0349;
    nlls::SVD<nlls::Matrix43d> svd(A);
    auto U = svd.U();
    auto S = svd.S();
    auto V = svd.V();
    REQUIRE(U(0, 0) == Approx(0.567329));
    REQUIRE(U(1, 2) == Approx(-0.726197));
    REQUIRE(U(3, 1) == Approx(-0.295112));
    REQUIRE(S[0] == Approx(6.10923));
    REQUIRE(S[1] == Approx(1.96659));
    REQUIRE(S[2] == Approx(1.39));
}