#pragma once

#include <type_traits>

namespace nlls {

template <typename T>
struct RemoveConstRef {
    using Type = typename std::remove_const<typename std::remove_reference<T>::type>::type;
};
template <typename T>
using RemoveConstRefType = typename RemoveConstRef<T>::Type;

template <bool B>
using EnableIfType = typename std::enable_if<B>::type;
template <bool B>
using EnableIfNotType = typename std::enable_if<!B>::type;

template <bool B, int N1, int N2>
struct IntSelector {
    static constexpr int N = B ? N1 : N2;
};

} // namespace nlls