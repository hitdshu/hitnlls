#include "node/node_pose2.h"
#include "catch/catch.hpp"

TEST_CASE( "Test node pose2", "[node]") {
    ::hitnlls::node::NodePose2 pos2;
    ::hitnlls::matrix::Vector3f update;
    update[0] = 0.7854;
    update[1] = 2;
    update[2] = 1;
    pos2.UpdatePlus(update);
    ::hitnlls::node::SE2 se2 = pos2.GetEstimate();
    REQUIRE(se2.ToMat33f()(0, 0) == Approx(0.707106));
    REQUIRE(se2.ToMat33f()(0, 1) == Approx(-0.707108));
    REQUIRE(se2.ToMat33f()(1, 0) == Approx(0.707108));
    REQUIRE(se2.ToMat33f()(1, 1) == Approx(0.707106));
}