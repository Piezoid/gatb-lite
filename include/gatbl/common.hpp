#ifndef COMMON_HPP
#define COMMON_HPP

/* GLIBC configuration */
#define _FILE_OFFSET_BITS 64

#ifndef PAGE_SIZE
/// Allow the user to use a global variable or a syscall for portability
#    define PAGE_SIZE (1UL << 12)
#endif

#include <cstdio>
#include <exception>

#define noinline_fun __attribute__((noinline))
#define forceinline_fun inline __attribute__((always_inline))
#define flatten_fun __attribute__((flatten))
#define pure_fun __attribute__((pure))
#define hot_fun __attribute__((hot))
#define cold_fun __attribute__((cold))
#define restrict __restrict__
#define likely(expr) __builtin_expect(!!(expr), 1)
#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define prefetchr(addr) __builtin_prefetch((addr), 0)
#define prefetchw(addr) __builtin_prefetch((addr), 1)

#ifndef __has_feature
#    define __has_feature(x) (0)
#endif

#ifdef __has_cpp_attribute
#    if __has_cpp_attribute(noreturn)
#        define noreturn_attr [[noreturn]]
#    endif
#    if __has_cpp_attribute(nodiscard)
#        define nodiscard_attr [[nodiscard]]
#    elif __has_cpp_attribute(gnu::warn_unused_result)
#        define nodiscard [[gnu::warn_unused_result]]
#    endif
#endif

#ifndef noreturn_attr
#    define noreturn_attr __attribute__((noreturn))
#endif

#ifndef nodiscard_attr
#    define nodiscard_attr __attribute__((warn_unused_result))
#endif

#ifdef NDEBUG
#    define DEBUG 0
#    define assert(expr) static_cast<void>(0)
#    define assume(expr) (likely((expr)) ? static_cast<void>(0) : __builtin_unreachable())
#else
#    define DEBUG 1

noreturn_attr noinline_fun inline void
              sourceloc_fail(const char* msg, const char* file, unsigned int line, const char* function) noexcept
{
    std::fprintf(stderr, "%s:%s:%u: %s\n", file, function, line, msg);
    std::fflush(stderr);
    std::terminate();
}

#    define assert(expr)                                                                                               \
        (likely((expr)) ? static_cast<void>(0)                                                                         \
                        : sourceloc_fail("Assertion '" #expr "' failed.", __FILE__, __LINE__, __PRETTY_FUNCTION__))

#    define assume assert

#endif

#endif // COMMON_HPP
