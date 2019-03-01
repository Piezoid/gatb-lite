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
    using element_type    = remove_reference_t<reference>;
    using repr_t          = I;

    iterator       begin() noexcept { return _begin; }
    const_iterator begin() const noexcept { return _begin; }
    sentinel       end() const noexcept { return _end; }

    template<typename _I = I, typename _S = S>
    auto size() const noexcept -> make_unsigned_t<decltype(std::declval<_S>() - std::declval<_I>())>
    {
        assume(this->end() >= this->begin(), "iterator ends before begin");
        return this->end() - this->begin();
    }
    bool empty() const noexcept { return !(this->begin() != this->end()); }

    template<typename _I = I> auto operator[](size_t i) const noexcept -> decltype(*(std::declval<_I>() + i))
    {
        assume(_begin + i < _end, "Out of bound access %lu >= %lu", i, size());
        return *(_begin + i);
    }

    template<typename _I = I> auto operator[](size_t i) noexcept -> decltype(*(std::declval<_I>() + i))
    {
        assume(_begin + i < _end, "Out of bound access %lu >= %lu", i, size());
        return *(_begin + i);
    }

    repr_t   _begin;
    sentinel _end;
};

template<typename I, typename S = I, typename CI = I>
auto
size(const iterator_pair<I, S, CI>& r) -> decltype(end(r) - begin(r))
{
    return end(r) - begin(r);
}

template<typename I, typename S, typename CI = I>
inline constexpr iterator_pair<I, S, CI>
make_range(I begin, S end)
{
    return {begin, end};
}

template<typename Container>
inline constexpr auto
make_range(Container& c)
  -> iterator_pair<decltype(begin(c)), decltype(end(c)), decltype(begin(static_cast<const Container&>(c)))>
{
    return {begin(c), end(c)};
}

} // namespace gatbl

#endif // INTERATOR_PAIR_HPP
