#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <cstdlib>

#include <exception>
#include <system_error>

// FIXME: separate assert.h
#include <sstream>
#include <iostream>

#include "gatbl/common.hpp"
#include "gatbl/utils/compatibility.hpp"

namespace gatbl { namespace sys {
// The oxymoric "noinline inline" means that we want a weak non inlineable function
noreturn_attr noinline_fun inline cold_fun void
throw_syserr(const char fmt[]...) // No variadic template, avoiding to generate too much code
{
    char msg[128];
    int  errcode = 0;
    std::swap(errcode, errno);
    va_list args;
    va_start(args, fmt);
    int size = vsnprintf(msg, sizeof(msg), fmt, args);
    va_end(args);
    if (size < 0 || size_t(size) >= sizeof(msg)) { std::terminate(); }
    throw std::system_error(errcode, std::generic_category(), msg);
}

template<typename T, typename... Args>
forceinline_fun hot_fun T
check_ret(T ret, T min, const char* what, Args&&... args)
{
    if (unlikely(ret < min)) { throw_syserr(what, std::forward<Args>(args)...); }
    return ret;
}

template<typename T, typename... Args>
forceinline_fun hot_fun T
check_ret(T ret, const char* what, Args&&... args)
{
    return check_ret(ret, T(0), what, std::forward<Args>(args)...);
}

template<typename T, typename... Args>
forceinline_fun hot_fun T*
                        check_ptr(T* ret, const char* what, Args&&... args)
{
    if (unlikely(ret == nullptr)) { throw_syserr(what, std::forward<Args>(args)...); }
    return ret;
}

} // namespace sys
} // namespace gatbl

#endif // EXCEPTIONS_H
