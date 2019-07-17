#include "node/so3.h"
#include "catch/catch.hpp"

TEST_CASE( "Test so3 initialization", "[node]") {
    hitnlls::node::SO3 so3;
    REQUIRE(so3.W() == 1);
    REQUIRE(so3.X() == 0);
    REQUIRE(so3.Y() == 0);
    REQUIRE(so3.Z() == 0);
}