#ifndef EMPTY_BASE_HPP
#define EMPTY_BASE_HPP

#include <utility>
#include "gatbl/utils/compatibility.hpp"

namespace gatbl { namespace utils {

/**
 * @brief Empty base optimization type
 *
 * In c++ an empty struct (eg. a stateless lambda) takes 1 byte. For many small objects this might
 * take considerable memory. A solution is to use the empty-base optimization where the empty type
 * is inherited by a structure containing a non empty data member. The compiler is allowed to size the struct
 * such as `sizeof(empty_base<T, Tag>) = sizeof(T)` iff `Tag` is an empty struct.
 */
template<typename T, typename Tag> class empty_base : private Tag
{
    using value_type = remove_reference_t<T>;

  public:
    template<typename _T = T, typename _Tag = Tag>
    empty_base(_T&& value, _Tag&& tag = Tag{})
      : Tag(std::forward<_Tag>(tag))
      , _value(std::forward<_T>(value))
    {}
    Tag&              tag() { return *this; }
    const Tag&        tag() const { return *this; }
    const value_type& value() const noexcept { return _value; }
    value_type&       value() noexcept { return _value; }

  private:
    T _value;
};

}}

#endif // EMPTY_BASE_HPP
