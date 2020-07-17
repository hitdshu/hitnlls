#pragma once

namespace nlls {
using std::abs;
using std::sqrt;
using std::pow;
using std::exp;
using std::log;
using std::sin;
using std::cos;
using std::asin;
using std::acos;
using std::atan2;
template <typename T>
class Jet {
public:
    using ThisType = Jet<T>;
    Jet() { val_ = 0; deriv_ = 0; }
    Jet(const T &v) { val_ = v; deriv_ = 0; }
    Jet(const T &v, const T &d) { val_ = v; deriv_ = d; }
    explicit operator T() const { return val_; }
    ThisType &operator+=(const ThisType &j) { val_ += j.val_; deriv_ += j.deriv_; return *this; }
    ThisType &operator-=(const ThisType &j) { val_ -= j.val_; deriv_ -= j.deriv_; return *this; }
    ThisType &operator*=(const ThisType &j) { val_ *= j.val_; deriv_ = deriv_*j.val_+j.deriv_*val_; return *this; }
    ThisType &operator/=(const ThisType &j) { val_ /= j.val_; deriv_ = (deriv_*j.val_-j.deriv_*val_)/(j.val_*j.val_); return *this; }
    ThisType &operator=(const T &v) { val_ = v; deriv_ = 0; }
    ThisType &operator+=(const T &v) { val_ += v; return *this; }
    ThisType &operator-=(const T &v) { val_ -= v; return *this; }
    ThisType &operator*=(const T &v) { val_ *= v; deriv_ *= v; return *this; }
    ThisType &operator/=(const T &v) { val_ /= v; deriv_ /= v; return *this; }
    bool operator>(const ThisType &j) const { return val_ > j.val_; }
    bool operator>=(const ThisType &j) const { return val_ >= j.val_; }
    bool operator<(const ThisType &j) const { return val_ < j.val_; }
    bool operator<=(const ThisType &j) const { return val_ <= j.val_; }
    bool operator==(const ThisType &j) const { return val_ == j.val_; }
    bool operator!=(const ThisType &j) const { return val_ != j.val_; }
    T &Val() { return val_; }
    const T &Val() const { return val_; }
    T &Deriv() { return deriv_; }
    const T &Deriv() const { return deriv_; }
private:
    T val_;
    T deriv_;
};
template <typename T>
Jet<T> operator+(const Jet<T> &j) {
    return j;
}
template <typename T>
Jet<T> operator-(const Jet<T> &j) {
    return Jet<T>(-j.Val(), -j.Deriv());
}
template <typename T>
Jet<T> operator+(const Jet<T> &j1, const Jet<T> &j2) {
    return Jet<T>(j1.Val()+j2.Val(), j1.Deriv()+j2.Deriv());
}
template <typename T>
Jet<T> operator+(const Jet<T> &j, const T &v) {
    return Jet<T>(j.Val()+v, j.Deriv());
}
template <typename T>
Jet<T> operator+(const T &v, const Jet<T> &j) {
    return j + v;
}
template <typename T>
Jet<T> operator-(const Jet<T> &j1, const Jet<T> &j2) {
    return Jet<T>(j1.Val()-j2.Val(), j1.Deriv()-j2.Deriv());
}
template <typename T>
Jet<T> operator-(const Jet<T> &j, const T &v) {
    return Jet<T>(j.Val()-v, j.Deriv());
}
template <typename T>
Jet<T> operator-(const T &v, const Jet<T> &j) {
    return Jet<T>(v-j.Val(), -j.Deriv());
}
template <typename T>
Jet<T> operator*(const Jet<T> &j1, const Jet<T> &j2) {
    return Jet<T>(j1.Val()*j2.Val(), j1.Deriv()*j2.Val()+j1.Val()*j2.Deriv());
}
template <typename T>
Jet<T> operator*(const T &v, const Jet<T> &j) {
    return Jet<T>(v*j.Val(), v*j.Deriv());
}
template <typename T>
Jet<T> operator*(const Jet<T> &j, const T &v) {
    return v * j;
}
template <typename T>
Jet<T> operator/(const T &v, const Jet<T> &j) {
    return Jet<T>(v/j.Val(), -v/j.Val()/j.Val()*j.Deriv());
}
template <typename T>
Jet<T> operator/(const Jet<T> &j, const T &v) {
    return j * (T(1.0) / v);
}
template <typename T>
Jet<T> operator/(const Jet<T> &j1, const Jet<T> &j2) {
    return Jet<T>(j1.Val()/j2.Val(), (j1.Deriv()*j2.Val()-j1.Val()*j2.Deriv())/(j2.Val()*j2.Val()));
}
template <typename T>
Jet<T> abs(const Jet<T> &j) {
    return Jet<T>(abs(j.Val()), (j.Val() > 0) ? j.Deriv() : -j.Deriv());
}
template <typename T>
Jet<T> sin(const Jet<T> &j) {
    return Jet<T>(sin(j.Val()), j.Deriv()*cos(j.Val()));
}
template <typename T>
Jet<T> cos(const Jet<T> &j) {
    return Jet<T>(cos(j.Val()), -j.Deriv()*sin(j.Val()));
}
template <typename T>
Jet<T> asin(const Jet<T> &j) {
    return Jet<T>(asin(j.Val()), j.Deriv()/sqrt(1-j.Val()*j.Val()));
}
template <typename T>
Jet<T> acos(const Jet<T> &j) {
    return Jet<T>(acos(j.Val()), -j.Deriv()/sqrt(1-j.Val()*j.Val()));
}
template <typename T>
Jet<T> exp(const Jet<T> &j) {
    return Jet<T>(exp(j.Val()), exp(j.Val())*j.Deriv());
}
template <typename T>
Jet<T> pow(const Jet<T> &j, int n) {
    return Jet<T>(pow(j.Val(),n), n*pow(j.Val(),n-1)*j.Deriv());
}
template <typename T>
Jet<T> log(const Jet<T> &j) {
    return Jet<T>(log(j.Val()), 1/j.Val()*j.Deriv());
}
template <typename T>
Jet<T> sqrt(const Jet<T> &j) {
    return Jet<T>(sqrt(j.Val()), 0.5/j.Val()*j.Deriv());
}
template <typename T>
Jet<T> atan2(const Jet<T> &j1, const Jet<T> &j2) {
    return Jet<T>(atan2(j1.Val(), j2.Val()), 1/(1+pow(j1.Val()/j2.Val(),2))*(j1.Deriv()*j2.Val()-j1.Val()*j2.Deriv())/(j1.Val()*j2.Val()));
}
template <class DestType, class SrcType>
void ConvertMatrix(DestType &dst, const SrcType &src) {
    typedef typename DestType::ElementType DestScalar;
    for (int i = 0; i < src.Rows(); ++i) {
        for (int j = 0; j < src.Cols(); ++j) {
            dst(i, j) = DestScalar(src(i, j));
        }
    }
}
template <typename T> 
std::ostream &operator<<(std::ostream &os, const Jet<T> &j) {
    os << "(" << j.Val() << "," << j.Deriv() << ")";
    return os;
}
using Jetd = Jet<double>;
using Jetf = Jet<float>;
} // namespace nlls