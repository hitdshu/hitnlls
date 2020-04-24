#pragma once

#include <map>
#include <algorithm>
#include <vector>

#include "nlls/trait.h"

namespace hitnlls {

template <typename T>
class SparseArray {
public:
    explicit SparseArray(int len = -1) { len_ = len; }
    SparseArray(const SparseArray &sa) { len_ = sa.len_; i2v_ = sa.i2v_; }

    bool CountIndex(int i) const { return i2v_.count(i); }
    int Length() const { return len_; }
    int NonzeroLen() const { return i2v_.size(); }
    void Delete(int i) { i2v_.erase(i); }
    int LengthAll() const {
        int total_len = 0;
        for (auto iter = i2v_.begin(); iter != i2v_.end(); ++iter) {
            total_len += LengthTraits<T>::Len(iter->second);
        }
        return total_len;
    }
    double Norm() const {
        double nsq = 0;
        for (auto iter = i2v_.begin(); iter != i2v_.end(); ++iter) {
            nsq += NormTraits<T>::NormSquare(iter->second);
        }
        return std::sqrt(nsq);
    }
    std::vector<int> GetIndices() const {
        std::vector<int> idxs;
        for (auto iter = i2v_.begin(); iter != i2v_.end(); ++iter) {
            idxs.push_back(iter->first);
        }
        return idxs;
    }

    T &operator()(int i) { return i2v_[i]; }
    T operator()(int i) const { auto iter = i2v_.find(i); if (iter == i2v_.end()) return T(); else return iter->second; }
    T &operator[](int i) { return i2v_[i]; }
    T operator[](int i) const { auto iter = i2v_.find(i); if (iter == i2v_.end()) return T(); else return iter->second; }
    SparseArray operator+(const SparseArray &other) const {
        SparseArray result(*this);
        result.len_ = std::max(len_, other.len_);
        std::vector<int> other_idxs = other.GetIndices();
        for (auto i : other_idxs) {
            if (i2v_.count(i)) {
                result[i] += other[i];
            } else {
                result[i] = other[i];
            }
        }
        return result;
    }
    SparseArray operator-(const SparseArray &other) const {
        SparseArray result(*this);
        result.len_ = std::max(len_, other.len_);
        std::vector<int> other_idxs = other.GetIndices();
        for (auto i : other_idxs) {
            if (i2v_.count(i)) {
                result[i] -= other[i];
            } else {
                result[i] = -other[i];
            }
        }
        return result;
    }
    SparseArray operator*(const T &val) const {
        SparseArray result(len_);
        std::vector<int> idxs = GetIndices();
        for (auto i : idxs) {
            result[i] = this->operator()(i) * val;
        }
        return result;
    }
    SparseArray &operator+=(const SparseArray &other) {
        std::vector<int> other_idxs = other.GetIndices();
        for (auto i : other_idxs) {
            if (i2v_.count(i)) {
                this->operator()(i) += other[i];
            } else {
                this->operator()(i) = other[i];
            }
        }
        return *this;
    }
    SparseArray &operator-=(const SparseArray &other) {
        std::vector<int> other_idxs = other.GetIndices();
        for (auto i : other_idxs) {
            if (i2v_.count(i)) {
                this->operator()(i) -= other[i];
            } else {
                this->operator()(i) = -other[i];
            }
        }
        return *this;
    }
    friend SparseArray operator*(const T &v, const SparseArray &sa) {
        SparseArray result(sa.len_);
        std::vector<int> idxs = sa.GetIndices();
        for (auto i : idxs) {
            result[i] = v * sa(i);
        }
        return result;
    }

private:
    std::map<int, T> i2v_;
    int len_;
};

using SparseArrayf = SparseArray<float>;
using SparseArrayd = SparseArray<double>;
using VectorSparseArray = SparseArray<matrix::VectorXd>;
using MatrixSparseArray = SparseArray<matrix::MatrixXd>;

} // namespace hitnlls