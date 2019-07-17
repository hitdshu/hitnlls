#include "node/so3.h"
#include "catch/catch.hpp"

TEST_CASE( "Test so3 initialization", "[node]") {
    hitnlls::node::SO3 so3;
    REQUIRE(so3.W() == Approx(1));
    REQUIRE(so3.X() == Approx(0));
    REQUIRE(so3.Y() == Approx(0));
    REQUIRE(so3.Z() == Approx(0));
}