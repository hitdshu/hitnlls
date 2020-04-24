#pragma once

namespace hitnlls {
namespace matrix {

struct Evaluator {
    template <class T, class S>
    HITNLLS_INLINE static void Assign(T &t, const S &s) {
        for (int i = 0; i < t.Rows(); ++i) {
            for (int j = 0; j < t.Cols(); ++j) {
                t.Cast()(i, j) = s.At(i, j);
            }
        }
    }

    template <class T, class S1, class S2>
    HITNLLS_INLINE static void Mult(T &t, const S1 &s1, const S2 &s2) {
        for (int i = 0; i < t.Rows(); ++i) {
            for (int j = 0; j < s2.Rows(); ++j) {
                for (int k = 0; k < t.Cols(); ++k) {
                    t.Cast()(i, k) += s1.At(i, j) * s2.At(j, k);
                }
            }
        }
    }
};

} // namespace matrix
} // namespace hitnlls