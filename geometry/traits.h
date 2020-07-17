#pragma once

#include "../matrix/dense.h"

namespace nlls {
template <typename T> class SO2;
template <typename T> class SO2Tangent;
template <typename T> class SE2;
template <typename T> class SE2Tangent;
template <typename T> class SO3;
template <typename T> class SO3Tangent;
template <typename T> class SE3;
template <typename T> class SE3Tangent;
namespace internal {

template <class GroupType> struct LieGroupTraits;
template <class GroupType> struct LieGroupTraits<const GroupType> : public LieGroupTraits<GroupType> {};
template <class GroupType> struct LieGroupTraits<GroupType &> : public LieGroupTraits<GroupType> {};
template <class TangentType> struct TangentTraits;
template <class TangentType> struct TangentTraits<const TangentType> : public TangentTraits<TangentType> {};
template <class TangentType> struct TangentTraits<TangentType &> : public TangentTraits<TangentType> {};

template <typename T> struct Constants { static constexpr float EPS = 1e-8; };
template <> struct Constants<double> { static constexpr double EPS = 1e-8; };
template <> struct Constants<float> { static constexpr float EPS = 1e-6; };
template <> struct Constants<Jetf> { static constexpr float EPS = 1e-6; };

} // namespace internal
} // namespace nlls