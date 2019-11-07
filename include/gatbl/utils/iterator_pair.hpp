#ifndef INTERATOR_PAIR_HPP
#define INTERATOR_PAIR_HPP

#include <string>

#include <iosfwd>
#include "gatbl/common.hpp"
#include "gatbl/utils/ranges.hpp"

namespace gatbl {

/// A Range made of an iterator and a sentinel, used as non owning proxy for a sequence of items
template<typename I, typename S = I> struct iterator_pair : public view_facade<iterator_pair<I, S>, I, I, S, S>
{
  private:
    using base = view_facade<iterator_pair<I, S>, I, I, S, S>;

  public:
    using typename base::size_type;

  protected:
    I _begin{};
    S _end{};

  public:
    iterator_pair() = default;
    iterator_pair(I beg, S end)
      : _begin(beg)
      , _end(end)
    {}

    template<typename R,
             typename = enable_if_t<!std::is_same<remove_reference_t<R>, iterator_pair>::value
                                    && std::is_convertible<iterator_t<R>, I>::value
                                    && std::is_convertible<sentinel_t<R>, S>::value>>
    constexpr iterator_pair(R&& range)
      : _begin(gatbl::begin(range))
      , _end(gatbl::end(range))
    {}

    template<typename _I, typename = decltype(concepts::value_require<I, S>(std::declval<_I>(), std::declval<_I>()))>
    CPP14_CONSTEXPR iterator_pair(const _I& it, size_type size)
      : _begin(it)
      , _end(it)
    {
        using std::advance;
        advance(_end, size);
    }

    iterator_pair(const iterator_pair&) = default;
    iterator_pair& operator=(const iterator_pair&) = default;

    I begin() const noexcept { return _begin; }
    S end() const noexcept { return _end; }
};

///// Specialization for pointer ranges (aka span)
// template<typename T> struct iterator_pair<T*, T*, const T*> {
//    using iterator        = T*;
//    using const_iterator  = const T*;
//    using sentinel        = T*;
//    using reference       = decltype(*std::declval<iterator>());
//    using const_reference = decltype(*std::declval<const_iterator>());
//    using value_type    = remove_reference_t<reference>;
//    using size_type       = make_unsigned_t<typename std::iterator_traits<iterator>::difference_type>;
//};

template<typename T> using span = iterator_pair<T*, T*>;

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
