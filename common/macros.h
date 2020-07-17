#pragma once

namespace nlls {

#define NLLS_WIN32 0x00
#define NLLS_LINUX 0x01

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    #define NLLS_PLATFORM NLLS_WIN32
    #define NLLS_INLINE __forceinline
#elif defined(__linux__)
    #define NLLS_PLATFORM NLLS_LINUX
    #define NLLS_INLINE __attribute__((always_inline))
#else
    #error "UNSUPPORTED PLATFORM"
#endif

#define NLLS_CRTP_REF \
    NLLS_INLINE const Derived &Cast() const { return static_cast<const Derived &>(*this); } \
    NLLS_INLINE Derived &Cast() { return static_cast<Derived &>(*this); }
#define NLLS_CRTP_DEC(T) \
    private: \
        NLLS_INLINE T() = default; \
        friend Derived;

#define NLLS_NONCOPYABLE(T) \
    T(const T &) = delete; \
    T(T &&) = delete; \
    T &operator=(const T &) = delete; \
    T &operator=(T &&) = delete;

} // namespace nlls