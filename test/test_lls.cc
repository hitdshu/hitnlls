#include "node/node_vec2.h"
#include "factor/factor_linear3_vec2.h"
#include "solver/gn_solver.h"
#include "graph/naive_graph.h"
#include "catch/catch.hpp"
#include <iostream>

TEST_CASE( "Test 1 factor/1 node lls", "[lls]") {
    ::hitnlls::node::NodeVec2 nv2;
    ::hitnlls::matrix::Vector2f x_est;
    x_est[0] = 0.7;
    x_est[1] = 0.3;
    nv2.SetEstimate(x_est);

    ::hitnlls::factor::FactorLinear3Vec2 fl3v2;
    ::hitnlls::matrix::Matrix32f matA;
    matA(0, 0) = 0.537667139546100;
    matA(0, 1) = 0.862173320368121;
    matA(1, 0) = 1.833885014595086;
    matA(1, 1) = 0.318765239858981;
    matA(2, 0) = -2.258846861003648;
    matA(2, 1) = -1.307688296305273;
    ::hitnlls::matrix::Vector3f b;
    b[0] = 0.700690668462573;
    b[1] = 1.882774391626265;
    b[2] = -2.155377560011862;
    fl3v2.SetA(matA);
    fl3v2.SetNode(0, &nv2);
    fl3v2.SetMeasurement(b);

    ::hitnlls::graph::NaiveGraph graph;
    graph.AddNode(&nv2);
    graph.AddFactor(&fl3v2);

    ::hitnlls::solver::GnSolver solver;
    solver.SetGraph(&graph);
    solver.Optimize();
    x_est = nv2.GetEstimate();

    REQUIRE(x_est[0] == Approx(0.964183));
    REQUIRE(x_est[1] == Approx(0.064255));
}