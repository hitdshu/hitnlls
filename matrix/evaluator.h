#pragma once

namespace nlls {
namespace internal {
struct Evaluator {
    template <class T, class S>
    NLLS_INLINE static void Assign(T &t, const S &s) {
        for (int i = 0; i < t.Rows(); ++i) {
            for (int j = 0; j < t.Cols(); ++j) {
                t.Cast()(i, j) = s.Cast()(i, j);
            }
        }
    }
    template <class T, class S1, class S2>
    NLLS_INLINE static void Mult(T &t, const S1 &s1, const S2 &s2) {
        for (int i = 0; i < t.Rows(); ++i) {
            for (int j = 0; j < s2.Rows(); ++j) {
                for (int k = 0; k < t.Cols(); ++k) {
                    t.Cast()(i, k) += s1.Cast()(i, j) * s2.Cast()(j, k);
                }
            }
        }
    }
};
} // namespace internal
} // namespace nlls