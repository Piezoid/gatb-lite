#ifndef PRECONDITION_CHECK_HPP
#define PRECONDITION_CHECK_HPP

#include "gatbl/common.hpp"

namespace gatbl { namespace utils {

template<bool debug = DEBUG> struct precondition_check;

template<> struct precondition_check<false>
{
    void set(bool) const {}
    void unchecked() const {}
    void checked() const {}
    void check() const {}
};

template<> struct precondition_check<true>
{
    void set(bool checked) const { _checked = checked; }
    void unchecked() const { _checked = false; }
    void checked() const { _checked = true; }
    void check() const { assume(_checked); }

  private:
    mutable bool _checked = false;
};

} // namespace utils
} // namespace gatbl

#endif // PRECONDITION_CHECK_HPP
