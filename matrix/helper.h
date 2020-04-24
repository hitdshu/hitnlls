#pragma once

namespace hitnlls {
namespace matrix {

struct Helper {
    template <typename ElementType>
    static bool CalcJacobRotCoeff(ElementType aii, ElementType aij, ElementType ajj, ElementType &cost, ElementType &sint) {
        if (std::abs(aij) < 1e-8) {
            cost = 1.0;
            sint = 0.0;
            return false;
        } else {
            ElementType tau = (ajj - aii) / (2 * aij);
            ElementType tmp = std::sqrt(1 + tau * tau);
            ElementType t1 = -tau + tmp;
            ElementType t2 = -tau - tmp;
            ElementType t = (std::abs(t1) > std::abs(t2)) ? t2 : t1;
            cost = 1 / std::sqrt(1 + t * t);
            sint = t * cost;
            return true;
        }
    }

    Helper() = delete;
};

} // namespace matrix
} // namespace hitnlls