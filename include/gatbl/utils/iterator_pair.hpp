#ifndef INTERATOR_PAIR_HPP
#define INTERATOR_PAIR_HPP

#include <string>

#include <iosfwd>
#include "gatbl/common.hpp"
#include "gatbl/utils/ranges.hpp"
#include "gatbl/utils/empty_base.hpp"

namespace gatbl {

/// A Range made of an iterator and a sentinel, used as non owning proxy for a sequence of items
template<typename I, typename S = I> struct iterator_pair : public view_facade<iterator_pair<I, S>, I, S>
{
  private:
    using it_traits = iterator_traits<I, S>;

  public:
    using iterator        = I;
    using sentinel        = S;
    using reference       = decltype(*std::declval<iterator>());
    using element_type    = remove_reference_t<reference>;
    using difference_type = typename it_traits::difference_type;
    using size_type       = typename it_traits::size_type;

  protected:
    iterator _begin;
    sentinel _end;

  public:
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

    template<typename _I,
             typename = decltype(concepts::value_require<I, S>(std::declval<_I>(), std::declval<_I>() + size_t(1)))>
    iterator_pair(const _I& it, size_t size)
      : _begin(it)
      , _end(it + size)
    {}

    iterator_pair(const iterator_pair&) = default;
    iterator_pair(iterator_pair&&)      = default;
    iterator_pair& operator=(const iterator_pair&) = default;
    iterator_pair& operator=(iterator_pair&&) = default;

    iterator begin() const noexcept { return _begin; }
    sentinel end() const noexcept { return _end; }
};

///// Specialization for pointer ranges (aka span)
// template<typename T> struct iterator_pair<T*, T*, const T*> {
//    using iterator        = T*;
//    using const_iterator  = const T*;
//    using sentinel        = T*;
//    using reference       = decltype(*std::declval<iterator>());
//    using const_reference = decltype(*std::declval<const_iterator>());
//    using element_type    = remove_reference_t<reference>;
//    using size_type       = make_unsigned_t<typename std::iterator_traits<iterator>::difference_type>;
//};

template<typename I, typename S, typename CI>
auto
size(const iterator_pair<I, S>& r) -> make_unsigned_t<decltype(distance(begin(r), end(r)))>
{
    return std::distance(begin(r), end(r));
}

template<typename T> using span = iterator_pair<T*, T*>;

template<typename T>
std::byte*
as_bytes(T* p)
{
    return reinterpret_cast<std::byte*>(p);
}

template<typename T>
const std::byte*
as_bytes(const T* p)
{
    return reinterpret_cast<const std::byte*>(p);
}

template<typename R>
auto
as_bytes(const R& r) -> decltype(span<const std::byte>(as_bytes(begin(r)), size(r) * sizeof(*begin(r))))
{
    return {as_bytes(begin(r)), size(r)};
}

template<typename R>
auto
as_bytes(R& r) -> decltype(span<std::byte>(as_bytes(begin(r)), size(r) * sizeof(*begin(r))))
{
    return {as_bytes(begin(r)), size(r) * sizeof(*begin(r))};
}

template<typename R> struct default_splitter
{
    using iterator = iterator_t<R>;
    using sentinel = sentinel_t<R>;

    explicit default_splitter(size_t size)
      : _size(size)
    {}
    default_splitter(const default_splitter&) = default;

    iterator split(iterator beg, sentinel end)
    {
        size_t dist = distance(beg, end);
        advance(beg, dist / 2);
        return dist > _size ? beg : end;
    }

    iterator_pair<iterator, sentinel> make_range(iterator beg, sentinel end) { return {beg, end}; }

    size_t _size = SIZE_MAX;
};

} // namespace gatbl

#endif // INTERATOR_PAIR_HPP
