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
ValueType DiagonalInverse(const ValueType &val) {
    return val.DiagonalInverse();
}
template <>
int DiagonalInverse<int>(const int &val) {
    return 1 / val;
}
template <>
float DiagonalInverse<float>(const float &val) {
    return 1.0 / val;
}
template <>
double DiagonalInverse<double>(const double &val) {
    return 1.0 / val;
}
template <typename ValueType>
struct DiagonalInverseTraits {
    typedef ValueType (*MDiaginverse)(const ValueType &);
    constexpr static MDiaginverse mdiaginv = DiagonalInverse<ValueType>;
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

template <typename ValueType>
double NormSquare(const ValueType &val) {
    double norm = val.Norm();
    return norm * norm;
}
template <>
double NormSquare<int>(const int &val) {
    return val * val;
}
template <>
double NormSquare<float>(const float &val) {
    return val * val;
}
template <>
double NormSquare<double>(const double &val) {
    return val * val;
}
template <typename ValueType>
struct NormSquareTraits {
    typedef double (*MNormsquare)(const ValueType &);
    constexpr static MNormsquare mnormsquare = NormSquare<ValueType>;
};

template <typename ValueType>
int Length(const ValueType &val) {
    return val.Rows();
}
template <>
int Length<int>(const int &val) {
    return 1;
}
template <>
int Length<float>(const float &val) {
    return 1;
}
template <>
int Length<double>(const double &val) {
    return 1;
}
template <typename ValueType>
struct LengthTraits {
    typedef int (*MLength)(const ValueType &);
    constexpr static MLength mlen = Length<ValueType>;
};
}

} // namespace matrix
} // namespace hitnlls