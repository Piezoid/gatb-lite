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
    using size_type       = make_unsigned_t<typename std::iterator_traits<iterator>::difference_type>;

    iterator_pair(I beg, S end)
      : _begin(beg)
      , _end(end)
    {}

    template<typename R,
             typename = enable_if_t<!std::is_same<remove_reference_t<R>, iterator_pair>::value
                                    && std::is_convertible<iterator_t<R>, iterator>::value
                                    && std::is_convertible<sentinel_t<R>, sentinel>::value>>
    iterator_pair(R& range)
      : _begin(gatbl::begin(range))
      , _end(gatbl::end(range))
    {}

    iterator_pair(const iterator_pair&) = default;
    iterator_pair(iterator_pair&&)      = default;
    iterator_pair& operator=(const iterator_pair&) = default;
    iterator_pair& operator=(iterator_pair&&) = default;

    iterator       begin() noexcept { return _begin; }
    const_iterator begin() const noexcept { return _begin; }
    sentinel       end() const noexcept { return _end; }

    template<typename _I = I, typename _S = S> auto size() const noexcept -> size_type
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

    iterator _begin;
    sentinel _end;
};

template<typename I, typename S, typename CI>
auto
size(const iterator_pair<I, S, CI>& r) -> make_unsigned_t<decltype(distance(begin(r), end(r)))>
{
    return std::distance(begin(r), end(r));
}

} // namespace gatbl

#endif // INTERATOR_PAIR_HPP
