#pragma once

namespace nlls {
namespace internal {
template <typename T, int N> 
class Blob {
public:
    NLLS_INLINE explicit Blob(int n) {}
    NLLS_INLINE explicit Blob() {}
    NLLS_INLINE explicit Blob(int n, const T &v) { Fill(v); }
    NLLS_INLINE explicit Blob(const T &v) { Fill(v); }
    NLLS_INLINE explicit Blob(const Blob &b) { data_ = b.data_; }
    NLLS_INLINE explicit Blob(Blob &&b) { data_ = std::move(b.data_); }
    NLLS_INLINE void Fill(const T &v) { data_.fill(v); }
    NLLS_INLINE void Swap(Blob &b) { data_.swap(b.data_); }
    NLLS_INLINE void Resize(int n) {}
    NLLS_INLINE Blob &operator=(const Blob &b) { data_ = b.data_; return *this; }
    NLLS_INLINE Blob &operator=(Blob &&b) { data_ = std::move(b.data_); return *this; }
    NLLS_INLINE T &operator[](int i) { return data_[i]; }
    NLLS_INLINE const T &operator[](int i) const { return data_[i]; }
    NLLS_INLINE T *Data() { return data_.data(); }
    NLLS_INLINE const T *Data() const { return data_.data(); }
    NLLS_INLINE int Size() const { return N; }
private:
    std::array<T, N> data_;
};
template <typename T>
class Blob<T, DYNAMIC> {
public:
    NLLS_INLINE explicit Blob(int n) : data_(n) {}
    NLLS_INLINE explicit Blob() {}
    NLLS_INLINE explicit Blob(int n, const T &v) : data_(n, v) {}
    NLLS_INLINE explicit Blob(const T &v) : data_(1, v) {}
    NLLS_INLINE explicit Blob(const Blob &b) { data_ = b.data_; }
    NLLS_INLINE explicit Blob(Blob &&b) { data_ = std::move(b.data_); }
    NLLS_INLINE void Fill(const T &v) { std::fill(data_.begin(), data_.end(), v); }
    NLLS_INLINE void Swap(Blob &b) { data_.swap(b.data_); }
    NLLS_INLINE void Resize(int n) { if (n != data_.size()) data_.resize(n); }
    NLLS_INLINE Blob &operator=(const Blob &b) { data_ = b.data_; return *this; }
    NLLS_INLINE Blob &operator=(Blob &&b) { data_ = std::move(b.data_); return *this; }
    NLLS_INLINE T &operator[](int i) { return data_[i]; }
    NLLS_INLINE const T &operator[](int i) const { return data_[i]; }
    NLLS_INLINE T *Data() { return data_.data(); }
    NLLS_INLINE const T *Data() const { return data_.data(); }
    NLLS_INLINE int Size() const { return data_.size(); }
private:
    std::vector<T> data_;
};
} // namespace internal
} // namespace nlls