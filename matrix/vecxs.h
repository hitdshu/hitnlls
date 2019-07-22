#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include "matrix/matrixx.h"
#include "matrix/mtraits.h"

namespace hitnlls {
namespace matrix {

template <typename ValueType>
class Vecxs {
public:
    Vecxs(int len = -1) { len_ = len; }
    Vecxs(const Vecxs<ValueType> &vec);

    inline ValueType &operator()(int idx);
    inline ValueType operator()(int idx) const;
    inline ValueType &operator[](int idx);
    inline ValueType operator[](int idx) const;
    Vecxs<ValueType> operator+(const Vecxs<ValueType> &vec) const;
    Vecxs<ValueType> operator-(const Vecxs<ValueType> &vec) const;
    ValueType operator*(const Vecxs<ValueType> &vec) const;
    Vecxs<ValueType> operator*(const ValueType &val) const;
    Vecxs<ValueType> &operator+=(const Vecxs<ValueType> &vec);
    Vecxs<ValueType> &operator+=(const ValueType &val);
    Vecxs<ValueType> &operator-=(const Vecxs<ValueType> &vec);
    Vecxs<ValueType> &operator-=(const ValueType &val);
    Vecxs<ValueType> &operator=(const Vecxs<ValueType> &vec);

    inline bool HasValue(int idx) const;
    inline void Delete(int idx);
    inline int Length() const { return len_; }
    inline int NzLen() const { return idx_value_.size(); }
    double Norm() const;
    void Print(const ::std::string &name = "") const;
    ::std::vector<int> GetIndices() const;
    ::std::vector<int> GetIndicesget(const int &idx) const;

    friend ::std::ostream &operator<<(::std::ostream &out, const Vecxs<ValueType> &vec) { vec.Print(); return out; }
    friend Vecxs<ValueType> operator*(const ValueType &val, const Vecxs<ValueType> &vec) {
        Vecxs<ValueType> result(vec.len_);
        ::std::vector<int> idxs = vec.GetIndices();
        for (size_t idx = 0; idx < idxs.size(); ++idx) {
            result[idx] = val * vec(idxs[idx]);   
        }
        return result;
    }

protected:


private:
    std::map<int, ValueType> idx_value_;
    int len_;
};

template <typename ValueType>
Vecxs<ValueType>::Vecxs(const Vecxs<ValueType> &vec) {
    len_ = vec.len_;
    idx_value_ = vec.idx_value_;
}

template <typename ValueType>
ValueType &Vecxs<ValueType>::operator()(int idx) {
    assert(idx < len_);
    return idx_value_[idx];
}

template <typename ValueType>
ValueType Vecxs<ValueType>::operator()(int idx) const {
    assert(idx < len_);
    auto iter = idx_value_.find(idx);
    if (iter != idx_value_.end()) {
        return iter->second;
    } else {
        ValueType ret_val;
        ret_val = 0;
        return ret_val;
    }
}

template <typename ValueType>
ValueType &Vecxs<ValueType>::operator[](int idx) {
    assert(idx < len_);
    return idx_value_[idx];
}

template <typename ValueType>
ValueType Vecxs<ValueType>::operator[](int idx) const {
    auto iter = idx_value_.find(idx);
    if (iter != idx_value_.end()) {
        return iter->second;
    } else {
        ValueType ret_val;
        ret_val = 0;
        return ret_val;
    }
}

template <typename ValueType>
bool Vecxs<ValueType>::HasValue(int idx) const {
    auto iter = idx_value_.find(idx);
    if (iter != idx_value_.end())
        return true;
    else 
        return false;
}

template <typename ValueType>
Vecxs<ValueType> Vecxs<ValueType>::operator+(const Vecxs<ValueType> &vec) const {
    Vecxs<ValueType> result(*this);
    result.len_ = ::std::max(len_, vec.len_);
    ::std::vector<int> idxs_vec = vec.GetIndices();
    for(size_t idx = 0; idx < idxs_vec.size(); ++idx) {
        if (idx_value_.count(idxs_vec[idx])) {
            result(idxs_vec[idx]) += vec(idxs_vec[idx]);
        } else {
            result(idxs_vec[idx]) = vec(idxs_vec[idx]);
        }
    }
    return result;
}

template <typename ValueType>
Vecxs<ValueType> Vecxs<ValueType>::operator-(const Vecxs<ValueType> &vec) const {
    Vecxs<ValueType> result(*this);
    result.len_ = ::std::max(len_, vec.len_);
    ::std::vector<int> idxs_vec = vec.GetIndices();
    for(size_t idx = 0; idx < idxs_vec.size(); ++idx) {
        if (idx_value_.count(idxs_vec[idx])) {
            result(idxs_vec[idx]) -= vec(idxs_vec[idx]);
        } else {
            result(idxs_vec[idx]) = - vec(idxs_vec[idx]);
        }
    }
    return result;
}

template <typename ValueType>
ValueType Vecxs<ValueType>::operator*(const Vecxs<ValueType> &vec) const {
    ValueType result(0);
    ::std::vector<int> idxs = GetIndices();
    for (size_t idx = 0; idx < idxs.size(); ++idx) {
        if (vec.idx_value_.count(idxs[idx])) {
            result += this->operator()(idxs[idx]) * vec[idxs[idx]];
        }
    }
    return result;
}

template <typename ValueType>
Vecxs<ValueType> Vecxs<ValueType>::operator*(const ValueType &val) const {
    Vecxs<ValueType> result(len_);
    ::std::vector<int> idxs = GetIndices();
    for (size_t idx = 0; idx < idxs.size(); ++idx) {
        result[idx] = this->operator()(idxs[idx]) * val;   
    }
    return result;
}

template <typename ValueType>
Vecxs<ValueType> &Vecxs<ValueType>::operator+=(const Vecxs<ValueType> &vec) {
    ::std::vector<int> idxs_vec = vec.GetIndices();
    for(size_t idx = 0; idx < idxs_vec.size(); ++idx) {
        if (idx_value_.count(idxs_vec[idx])) {
            this->operator()(idxs_vec[idx]) += vec[idxs_vec[idx]];
        } else {
            this->operator()(idxs_vec[idx]) = vec[idxs_vec[idx]];
        }
    }
    return *this;
}

template <typename ValueType>
Vecxs<ValueType> &Vecxs<ValueType>::operator+=(const ValueType &val) {
    ::std::vector<int> idxs = GetIndices();
    for (size_t idx = 0; idx < idxs.size(); ++idx) {
        this->operator()(idxs[idx]) += val;
    }
    return *this;
}

template <typename ValueType>
Vecxs<ValueType> &Vecxs<ValueType>::operator-=(const Vecxs<ValueType> &vec) {
    ::std::vector<int> idxs_vec = vec.GetIndices();
    for(size_t idx = 0; idx < idxs_vec.size(); ++idx) {
        if (idx_value_.count(idxs_vec[idx])) {
            this->operator()(idxs_vec[idx]) -= vec[idxs_vec[idx]];
        } else {
            this->operator()(idxs_vec[idx]) = - vec[idxs_vec[idx]];
        }
    }
    return *this;
}

template <typename ValueType>
Vecxs<ValueType> &Vecxs<ValueType>::operator-=(const ValueType &val) {
    ::std::vector<int> idxs = GetIndices();
    for(size_t idx = 0; idx < idxs.size(); ++idx) {
        if (idx_value_.count(idxs[idx])) {
            this->operator()(idxs[idx]) -= val;
        }
    }
    return *this;
}

template <typename ValueType>
Vecxs<ValueType> &Vecxs<ValueType>::operator=(const Vecxs<ValueType> &vec) {
    len_ = vec.len_;
    idx_value_ = vec.idx_value_;
    return *this;
}

template <typename ValueType>
void Vecxs<ValueType>::Delete(int idx) {
    idx_value_.erase(idx);
}

template <typename ValueType>
void Vecxs<ValueType>::Print(const ::std::string &name) const {
    if (0 == idx_value_.size()) {
        ::std::cout << ::std::endl;
        return;
    }
    ::std::cout << "Sparse vector " << name << " with length " << len_ << " "; 
    for (auto iter = idx_value_.begin(); iter != idx_value_.end(); ++iter) {
        ::std::cout << "(" << iter->first << ", " << iter->second << ") ";
    }
    ::std::cout << ::std::endl;
}

template <typename ValueType>
::std::vector<int> Vecxs<ValueType>::GetIndices() const {
    ::std::vector<int> idxs;
    for (auto iter = idx_value_.begin(); iter != idx_value_.end(); ++iter) {
        idxs.push_back(iter->first);
    }
    return idxs;
}

template <typename ValueType>
::std::vector<int> Vecxs<ValueType>::GetIndicesget(const int &idx) const {
    ::std::vector<int> idxs;
    for (auto iter = idx_value_.begin(); iter != idx_value_.end(); ++iter) {
        if (iter->first >= idx) {
            idxs.push_back(iter->first);
        }
    }
    return idxs;
}

template <typename ValueType>
double Vecxs<ValueType>::Norm() const {
    double norm_square = 0;
    for (auto iter = idx_value_.begin(); iter != idx_value_.end(); ++iter) {
        norm_square += NormSquareTraits<ValueType>::mnormsquare(iter->second);
    }
    return sqrt(norm_square);
}

using Vecis = Vecxs<int>;
using Vecfs = Vecxs<float>;
using Vecds = Vecxs<double>;

using Vecxsxf = Vecxs<Matrixxf>;
using Vecxsxd = Vecxs<Matrixxd>;
using Vecxsxi = Vecxs<Matrixxi>;

} // namespace matrix
} // namespace hitnlls