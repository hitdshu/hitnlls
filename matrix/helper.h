#pragma once

namespace nlls {
namespace internal {
template <typename ElementType>
bool CalcJacobRotCoeff(ElementType aii, ElementType aij, ElementType ajj, ElementType &cost, ElementType &sint) {
    if (nlls::abs(aij) < ElementType(1e-8)) {
        cost = ElementType(1.0);
        sint = ElementType(0.0);
        return false;
    } else {
        ElementType tau = (ajj - aii) / (ElementType(2) * aij);
        ElementType tmp = nlls::sqrt(ElementType(1.0) + tau * tau);
        ElementType t1 = -tau + tmp;
        ElementType t2 = -tau - tmp;
        ElementType t = (nlls::abs(t1) > nlls::abs(t2)) ? t2 : t1;
        cost = ElementType(1.0) / nlls::sqrt(ElementType(1.0) + t * t);
        sint = t * cost;
        return true;
    }
}
template <typename T>
T RandomNormal() {
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution<double> d{0, 1};
    return T(d(gen));
}
} // namespace internal
} // namespace nlls