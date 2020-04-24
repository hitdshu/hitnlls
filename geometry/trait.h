#pragma once

#include "matrix/dense.h"

namespace hitnlls {
namespace geometry {

template <class GroupType> struct LieGroupTraits;
template <class TangentType> struct TangentTraits;

template <typename T> struct Constants {
    static constexpr double EPS = 1e-8;
};
template <>
struct Constants<double> {
    static constexpr double EPS = 1e-8;
};
template <>
struct Constants<float> {
    static constexpr float EPS = 1e-8;
};

} // namespace geometry
} // namespace hitnlls