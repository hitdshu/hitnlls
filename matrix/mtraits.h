#pragma once

#include <cmath>

namespace hitnlls {
namespace matrix {

namespace {
template <typename ValueType>
ValueType CholeskyLLT(const ValueType &val) {
    return val.CholeskyLLT();
}

template <>
int CholeskyLLT<int>(const int &val) {
    return sqrt(val);
}

template <>
float CholeskyLLT<float>(const float &val) {
    return sqrt(val);
}

template <>
double CholeskyLLT<double>(const double &val) {
    return sqrt(val);
}

template <typename ValueType>
struct MsqrtTraits {
    typedef ValueType (*Msqrt)(const ValueType &);
    constexpr static Msqrt msqrt = CholeskyLLT<ValueType>;
};

template <typename ValueType>
ValueType Transpose(const ValueType &val) {
    return val.Transpose();
}

template <>
int Transpose<int>(const int &val) {
    return val;
}

template <>
float Transpose<float>(const float &val) {
    return val;
}

template <>
double Transpose<double>(const double &val) {
    return val;
}

template <typename ValueType>
struct TransposeTraits {
    typedef ValueType (*Mtranspose)(const ValueType &);
    constexpr static Mtranspose mtranspose = Transpose<ValueType>;
};

template <typename ValueType>
ValueType Inverse(const ValueType &val) {
    return val.Inverse();
}

template <>
int Inverse<int>(const int &val) {
    return 1 / val;
}

template <>
float Inverse<float>(const float &val) {
    return 1.0 / val;
}

template <>
double Inverse<double>(const double &val) {
    return 1.0 / val;
}

template <typename ValueType>
struct InverseTraits {
    typedef ValueType (*Minverse)(const ValueType &);
    constexpr static Minverse minv = Inverse<ValueType>;
};

template <typename ValueType>
void SetLowtri(ValueType &val) {
    val.SetLowtri();
}

template <>
void SetLowtri<int>(int &val) {
    return;
}

template <>
void SetLowtri<float>(float &val) {
    return;
}

template <>
void SetLowtri<double>(double &val) {
    return;
}

template <typename ValueType>
struct SetLowtriTraits {
    typedef void (*MSetLowtri)(ValueType &);
    constexpr static MSetLowtri msetlt = SetLowtri<ValueType>;
};

template <typename ValueType>
ValueType Identity(const ValueType &val) {
    return ValueType::Identity(val.Rows());
}

template <>
int Identity<int>(const int &val) {
    return 1;
}

template <>
float Identity<float>(const float &val) {
    return 1.0;
}

template <>
double Identity<double>(const double &val) {
    return 1.0;
}

template <typename ValueType>
struct IdentityTraits {
    typedef ValueType (*MIdentity)(const ValueType &);
    constexpr static MIdentity midentity = Identity<ValueType>;
};
}

} // namespace matrix
} // namespace hitnlls