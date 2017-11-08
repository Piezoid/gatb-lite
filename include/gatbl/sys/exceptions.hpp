#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <cstdlib>

#include <exception>
#include <system_error>

#include "gatbl/common.hpp"

namespace gatbl { namespace sys {
// The oxymoric "noinline inline" means that we want a weak non inlineable function
noreturn_attr noinline_fun inline cold_fun void
              throw_syserr(const char* what)
{
    int errcode = errno;
    errno       = 0;
    throw std::system_error(errcode, std::generic_category(), what);
}

template<typename T, typename Tbounded>
forceinline_fun hot_fun Tbounded
                        check_ret(T ret, Tbounded min, const char* what)
{
    if (unlikely(ret < T(min))) { throw_syserr(what); }
    return Tbounded(ret);
}

template<typename T>
forceinline_fun hot_fun std::make_unsigned_t<T>
                        check_ret(T ret, const char* what)
{
    return check_ret(ret, std::make_unsigned_t<T>(0), what);
}

template<typename T>
forceinline_fun hot_fun T*
                        check_ptr(T* ret, const char* what)
{
    if (unlikely(ret == nullptr)) { throw_syserr(what); }
    return ret;
}

} // namespace sys
} // namespace gatbl

#endif // EXCEPTIONS_H
