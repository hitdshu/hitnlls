#pragma once

#include "matrix/dense.h"

namespace hitnlls {

template <typename T>
struct LengthTraits {
    static inline int Len(const T &v) { return v.Rows(); }
};

template <>
struct LengthTraits<int> {
    static inline int Len(const int &v) { return 1; }
};

template <>
struct LengthTraits<float> {
    static inline int Len(const float &v) { return 1; }
};

template <>
struct LengthTraits<double> {
    static inline int Len(const double &v) { return 1; }
};

template <typename T>
struct NormTraits {
    static double NormSquare(const T &v) { double n = v.Norm(); return n*n; }
};

template <>
struct NormTraits<int> {
    static double NormSquare(const int &v) { return v*v; }
};

template <>
struct NormTraits<float> {
    static double NormSquare(const float &v) { return v*v; }
};

template <>
struct NormTraits<double> {
    static double NormSquare(const double &v) { return v*v; }
};

template <typename T>
struct IdentityTraits {
    static T Identity(const T &v) { T result = v; result.SetIdentity(); return result; }
};

template <>
struct IdentityTraits<int> {
    static int Identity(const int &v) { return 1; }
};

template <>
struct IdentityTraits<float> {
    static float Identity(const float &v) { return 1.0; }
};

template <>
struct IdentityTraits<double> {
    static double Identity(const double &v) { return 1.0; }
};

template <typename T>
struct SqrtTraits {
    static T Sqrt(const T &v) { return v.Sqrt(); }
};

template <>
struct SqrtTraits<int> {
    static int Sqrt(const int &v) { return std::sqrt(v); }
};

template <>
struct SqrtTraits<float> {
    static float Sqrt(const float &v) { return std::sqrt(v); }
};

template <>
struct SqrtTraits<double> {
    static double Sqrt(const double &v) { return std::sqrt(v); }
};

template <typename T>
struct InverseTraits {
    static T Inverse(const T &v) { return v.Inverse(); }
};

template <>
struct InverseTraits<int> {
    static int Inverse(const int &v) { return 1 / v; }
};

template <>
struct InverseTraits<float> {
    static float Inverse(const float &v) { return 1 / v; }
};

template <>
struct InverseTraits<double> {
    static double Inverse(const double &v) { return 1 / v; }
};

template <typename T>
struct TransposeTraits {
    static T Transpose(const T &v) { return v.Transpose(); }
};

template <>
struct TransposeTraits<int> {
    static int Transpose(const int &v) { return v; }
};

template <>
struct TransposeTraits<float> {
    static float Transpose(const float &v) { return v; }
};

template <>
struct TransposeTraits<double> {
    static double Transpose(const double &v) { return v; }
};

} // namespace hitnlls