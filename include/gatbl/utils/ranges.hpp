#ifndef RANGES_HPP
#define RANGES_HPP

#include <cstring>
#include <utility>
#include <memory>
#include <iterator>
#include <functional>
#include <algorithm>
#include <ostream>
#include "gatbl/utils/concepts.hpp"
#include "gatbl/utils/compatibility.hpp"

namespace gatbl {

/// Customization point for random acces
template<typename Range>
auto
at(Range& r, size_t pos) -> decltype(concepts::value_require(*(begin(r) + pos), size(r)))
{
    assume(pos <= size(r), "out of bound i=%lu < size()=%lu", pos, size(r));
    return *(begin(r) + pos);
}

// Those are not customization points but we will overload them with range-like arguments
using std::equal;
using std::find;
using std::lexicographical_compare;

namespace concepts {

template<typename R>
auto
Range(R&& r) -> decltype(begin(r) != end(r));

template<typename R>
auto
FiniteRange(R&& r) -> decltype(value_require(begin(r) != end(r), size(r)));

template<typename R, typename T>
auto
RangeOf(R&& r) -> decltype(type_require<T>(begin(r) != end(r), static_cast<T&>(*begin(r))));

template<typename R, typename T>
auto
FiniteRangeOf(R&& r) -> decltype(type_require<T>(begin(r) != end(r), size(r), static_cast<T&>(*begin(r))));

template<typename... T>
auto
Ranges(T&&... rs) -> decltype(value_require(Range(rs)...));

template<typename R>
auto
CStringRange(R& r) -> decltype(value_require<const char*>(r.begin()), size(r));

template<typename X, typename Y, typename Zipper>
auto
ZippableRanges(X& x, Y& y, Zipper&& binop = Zipper{})
  -> decltype(value_require(Ranges(x, y), binop(*begin(x), *begin(y))));

}

template<typename To, typename From>
inline constexpr auto
convert(From&& x) -> decltype(concepts::value_require<To>(std::forward<From>(x)))
{
    return std::forward<From>(x);
}

template<typename Range> using iterator_t         = decltype(begin(std::declval<Range&>()));
template<typename Range> using sentinel_t         = decltype(end(std::declval<Range&>()));
template<typename Range> using reverse_iterator_t = decltype(rbegin(std::declval<Range&>()));
template<typename Range> using reverse_sentinel_t = decltype(rend(std::declval<Range&>()));
template<typename Range> using reference_t        = decltype(*begin(std::declval<Range&>()));
template<typename Range> using const_reference_t  = decltype(*begin(std::declval<const Range&>()));
template<typename Range> using value_t            = remove_reference_t<reference_t<Range>>;

// [string.view.io], Inserters and extractors
template<typename R, typename Traits, typename = decltype(concepts::CStringRange<const R>)>
inline std::basic_ostream<char, Traits>&
operator<<(std::basic_ostream<char, Traits>& os, const R& s)
{
    return os.write(begin(s), size(s));
}

/// Find specialization for char like types
template<typename T, typename = enable_if_t<sizeof(T) == 1>>
inline T*
find(T* first, T* last, const remove_const_t<T>& v)
{
    T* p = static_cast<T*>(memchr(first, reinterpret_cast<const int&>(v), last - first));
    return likely(p != nullptr) ? p : last;
}

template<typename Range, typename T>
inline auto
find(const Range& r, const T& v) -> decltype(find(begin(r), end(r), v))
{
    return find(begin(r), end(r), v);
}

template<typename I, typename S, typename T>
inline auto
contains(I f, S l, const T& v) -> decltype(find(f, l, v) != l)
{
    return find(f, l, v) != l;
}

template<typename Range, typename T>
inline auto
contains(const Range& r, const T& v) -> decltype(find(r, v) != end(r))
{
    return find(r, v) != end(r);
}

template<typename X, typename Y, typename C, typename = decltype(concepts::Ranges<X, Y>)>
inline constexpr bool
equal(const X& x, const Y& y) noexcept
{
    return equal(begin(x), end(y), begin(y), end(y));
}

template<typename X, typename Y, typename = decltype(concepts::ZippableRanges<X, Y, std::less<X>>)>
inline constexpr bool
lexicographical_compare(const X& x, const Y& y)
{
    return lexicographical_compare(begin(x), end(x), begin(y), end(y));
}

template<typename X, typename Y, typename C, typename = decltype(concepts::ZippableRanges<X, Y, C>)>
inline constexpr bool
lexicographical_compare(const X& x, const Y& y, C&& comp)
{
    return lexicographical_compare(begin(x), end(x), begin(y), end(y), std::forward<C>(comp));
}

// template<typename X,typename Y, typename=decltype(concepts::ZippableRanges<X, Y, std::equal_to<>>)>
// inline constexpr bool operator==(const X& x, const Y& y) noexcept { return equal(x, y); }

// template<typename X,typename Y, typename=decltype(concepts::ZippableRanges<X, Y, std::equal_to<>>)>
// inline constexpr bool operator!=(const X& x, const Y& y) noexcept { return !equal(x, y); }

// template<typename X,typename Y, typename=decltype(concepts::ZippableRanges<X, Y, std::less<>>)>
// inline constexpr bool operator<(const X& x, const Y& y) noexcept { return lexicographical_compare(x, y); }

// template<typename X,typename Y, typename=decltype(concepts::ZippableRanges<X, Y, std::less_equal<>>)>
// inline constexpr bool operator<=(const X& x, const Y& y) noexcept { return lexicographical_compare(x, y,
// std::less_equal{}); }

// template<typename X,typename Y, typename=decltype(concepts::ZippableRanges<X, Y, std::greater<>>)>
// inline constexpr bool operator>(const X& x, const Y& y) noexcept { return lexicographical_compare(x, y,
// std::greater{}); }

// template<typename X,typename Y, typename=decltype(concepts::ZippableRanges<X, Y, std::greater_equal<>>)>
// inline constexpr bool operator>=(const X& x, const Y& y) noexcept { return lexicographical_compare(x, y,
// std::greater_equal{}); }

}

#endif // RANGES_HPP
