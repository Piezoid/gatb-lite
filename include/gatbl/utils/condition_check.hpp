#ifndef condition_check_HPP
#define condition_check_HPP

#include "gatbl/common.hpp"

namespace gatbl { namespace utils {

#ifdef NDEBUG

struct condition_check
{
    void set(bool) const {}
    void unchecked() const {}
    void checked() const {}
    void check() const {}
};

#else

struct condition_check
{
    void set(bool checked) const { _checked = checked; }
    void unchecked() const { _checked = false; }
    void checked() const { _checked = true; }
    void check(const char file[] = __builtin_FILE(),
               const char fun[]  = __builtin_FUNCTION(),
               unsigned   line   = __builtin_LINE()) const
    {
        if (unlikely(!_checked)) { abort_message("%s:%s:%u: condition check failed", file, fun, line); }
    }

  private:
    mutable bool _checked = false;
};

#endif

} // namespace utils
} // namespace gatbl

#endif // condition_check_HPP
