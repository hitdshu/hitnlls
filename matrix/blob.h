#pragma once

namespace hitnlls {
namespace matrix {

template <typename T, int N> 
class Blob {
public:
    HITNLLS_INLINE explicit Blob(int n) {}
    HITNLLS_INLINE explicit Blob() {}
    HITNLLS_INLINE explicit Blob(int n, const T &v) { Fill(v); }
    HITNLLS_INLINE explicit Blob(const T &v) { Fill(v); }
    HITNLLS_INLINE explicit Blob(const Blob &b) { data_ = b.data_; }

    HITNLLS_INLINE void Fill(const T &v) { data_.fill(v); }
    HITNLLS_INLINE void Swap(Blob &b) { data_.swap(b.data_); }
    HITNLLS_INLINE void Resize(int n) {}

    HITNLLS_INLINE Blob &operator=(const Blob &b) { data_ = b.data_; return *this; }
    HITNLLS_INLINE T &operator[](int i) { return data_[i]; }
    HITNLLS_INLINE const T &operator[](int i) const { return data_[i]; }

    HITNLLS_INLINE T *Data() { return data_.data(); }
    HITNLLS_INLINE const T *Data() const { return data_.data(); }

    HITNLLS_INLINE int Size() const { return N; }

private:
    std::array<T, N> data_;
};

template <typename T>
class Blob<T, DYNAMIC> {
public:
    HITNLLS_INLINE explicit Blob(int n) : data_(n) {}
    HITNLLS_INLINE explicit Blob() {}
    HITNLLS_INLINE explicit Blob(int n, const T &v) : data_(n, v) {}
    HITNLLS_INLINE explicit Blob(const T &v) : data_(1, v) {}
    HITNLLS_INLINE explicit Blob(const Blob &b) { data_ = b.data_; }
    HITNLLS_INLINE explicit Blob(Blob &&b) { data_ = std::move(b.data_); }

    HITNLLS_INLINE void Fill(const T &v) { std::fill(data_.begin(), data_.end(), v); }
    HITNLLS_INLINE void Swap(Blob &b) { data_.swap(b.data_); }
    HITNLLS_INLINE void Resize(int n) { if (n != data_.size()) data_.resize(n); }

    HITNLLS_INLINE Blob &operator=(const Blob &b) { data_ = b.data_; return *this; }
    HITNLLS_INLINE Blob &operator=(Blob &&b) { data_ = std::move(b.data_); return *this; }
    HITNLLS_INLINE T &operator[](int i) { return data_[i]; }
    HITNLLS_INLINE const T &operator[](int i) const { return data_[i]; }

    HITNLLS_INLINE T *Data() { return data_.data(); }
    HITNLLS_INLINE const T *Data() const { return data_.data(); }

    HITNLLS_INLINE int Size() const { return data_.size(); }

private:
    std::vector<T> data_;
};

} // namespace matrix
} // namespace hitnlls