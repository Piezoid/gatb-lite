#ifndef INTERATOR_PAIR_HPP
#define INTERATOR_PAIR_HPP

#include "gatbl/common.hpp"
#include "gatbl/utils/ranges.hpp"
#include "gatbl/utils/empty_base.hpp"

namespace gatbl {

/// A Range made of an iterator and a sentinel
template<typename I, typename S = I, typename CI = I> struct iterator_pair
{
    using iterator        = I;
    using const_iterator  = CI;
    using sentinel        = S;
    using reference       = decltype(*std::declval<iterator>());
    using const_reference = decltype(*std::declval<const_iterator>());
    using element_type    = std::remove_reference_t<reference>;
    using repr_t          = I;

    iterator       begin() noexcept { return _begin; }
    const_iterator begin() const noexcept { return _begin; }
    sentinel       end() const noexcept { return _end; }

    auto size() const noexcept -> decltype(this->end() - this->begin())
    {
        assume(this->end() >= this->begin());
        return this->end() - this->begin();
    }
    bool empty() const noexcept { return !(this->begin() != this->end()); }

    const_reference operator[](size_t i) const noexcept
    {
        assume(_begin + i < _end);
        return *(_begin + i);
    }

    reference operator[](size_t i) noexcept
    {
        assume(_begin + i < _end);
        return *(_begin + i);
    }

    repr_t   _begin = {};
    sentinel _end   = {};
};

template<typename Container>
inline constexpr auto
make_range(Container& c)
  -> iterator_pair<decltype(begin(c)), decltype(end(c)), decltype(begin(static_cast<const Container&>(c)))>
{
    return {begin(c), end(c)};
}

} // namespace gatbl

#endif // INTERATOR_PAIR_HPP
