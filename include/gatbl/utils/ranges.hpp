#ifndef GATBL_RANGES_HPP
#define GATBL_RANGES_HPP

#include <cstring>
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
Range(R&& r) -> decltype(type_require<R>(begin(r) != end(r)));

template<typename R>
auto
FiniteRange(R&& r) -> decltype(type_require<R>(begin(r) != end(r), size(r)));

template<typename T, typename R>
auto
RangeOf(R&& r) -> decltype(type_require<R>(begin(r) != end(r), value_require<T&>(*begin(r))));

template<typename T, typename R>
auto
FiniteRangeOf(R&& r) -> decltype(type_require<T>(begin(r) != end(r), size(r), static_cast<T&>(*begin(r))));

template<typename... T>
auto
Ranges(T&&... rs) -> decltype(value_require(Range(rs)...));

template<typename R>
auto
CStringRange(R& r) -> decltype(value_require<const char*>(r.begin()), size(r));

template<typename Zipper, typename X, typename Y>
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

template<typename Iterator,
         typename Sentinel = Iterator,
         typename Category = typename std::iterator_traits<Iterator>::iterator_category>
struct iterator_traits : public std::iterator_traits<Iterator>
{
  protected:
    using base   = std::iterator_traits<Iterator>;
    using traits = iterator_traits<Iterator, Sentinel>; // get the most specialized version
    using IR     = Iterator const&;
    using SR     = Sentinel const&;

  public:
    using iterator = Iterator;
    using sentinel = Sentinel;
    using typename base::difference_type;
    using typename base::reference;
    using typename base::value_type;
    using size_type = make_unsigned_t<difference_type>;

    static constexpr bool is_output        = !std::is_const<value_type>::value;
    static constexpr bool is_bidirectional = false;
    static constexpr bool is_random_access = false;

    static bool      bound_check(IR it, IR, SR last) { return it != last; }
    static reference deref_check(IR it, IR first, SR last)
    {
        assume(traits::bound_check(it, first, last), "Out of bound");
        return *it;
    }
    static reference front(IR first, SR last) { return traits::deref_check(first, first, last); }
};

template<typename I, typename S>
struct iterator_traits<I, S, std::bidirectional_iterator_tag> : public iterator_traits<I, S, std::forward_iterator_tag>
{
  protected:
    using base = iterator_traits<I, S, std::forward_iterator_tag>;
    using typename base::IR;
    using typename base::SR;
    using typename base::traits;

  public:
    using typename base::reference;
    static constexpr bool is_bidirectional = true;

    static reference back(IR first, SR last)
    {
        auto it = last;
        return traits::deref_check(--it, first, last);
    }
};

template<typename I, typename S>
struct iterator_traits<I, S, std::random_access_iterator_tag>
  : public iterator_traits<I, S, std::bidirectional_iterator_tag>
{
  protected:
    using base = iterator_traits<I, S, std::bidirectional_iterator_tag>;
    using typename base::IR;
    using typename base::SR;
    using typename base::traits;

  public:
    using typename base::difference_type;
    using typename base::reference;
    using typename base::size_type;
    static constexpr bool is_random_access = true;

    static bool bound_check(IR it, IR first, SR last)
    {
        assume(first <= last, "iterator ends swapped");
        return first <= it && it < last;
    }
    static size_type size(IR first, SR last)
    {
        difference_type delta = last - first;
        assume(delta >= 0, "iterator ends swapped");
        return static_cast<size_type>(delta);
    }
    static reference at(size_type idx, IR first, SR last) { return traits::deref_check(first + idx, first, last); }
};

// using plopi = decltype(iterator_traits<int*, int*>::at);

template<typename Derived> struct iterator_facade
{
  protected:
    template<typename D = Derived, typename = enable_if_t<std::is_same<D, Derived>::value>> D& derived()
    {
        return *static_cast<Derived*>(this);
    }
    template<typename D = Derived, typename = enable_if_t<std::is_same<D, Derived>::value>> const D& derived() const
    {
        return *static_cast<const Derived*>(this);
    }

  private:
    template<typename D> auto operator*() -> decltype(derived<D>().dereference())
    {
        return this->derived().dereference();
    }

    template<typename D> auto operator++() -> decltype(concepts::type_require<Derived&>(derived<D>().increment()))
    {
        this->derived().increment();
        return this->derived();
    }

    template<typename D> auto operator--() -> decltype(concepts::type_require<Derived&>(derived<D>().decrement()))
    {
        this->derived().decrement();
        return this->derived();
    }

    template<typename D>
    auto operator++(int) -> decltype(concepts::type_require<Derived>(derived<D>().increment(), *this = *this))
    {
        Derived tmp(this->derived());
        this->derived().increment();
        return tmp;
    }

    template<typename D>
    auto operator--(int) -> decltype(concepts::type_require<Derived>(derived<D>().decrement(), *this = *this))
    {
        Derived tmp(this->derived());
        this->derived().decrement();
        return tmp;
    }

    template<typename D, typename d_t>
    auto operator+=(d_t n) -> decltype(concepts::type_require<Derived>(derived<D>().advance(n)))
    {
        this->derived().advance(n);
        return this->derived();
    }

    template<typename D, typename d_t>
    auto operator-=(d_t n) -> decltype(concepts::type_require<Derived>(derived<D>().advance(-n)))
    {
        this->derived().advance(-n);
        return this->derived();
    }
};

template<typename I, typename S, typename T>
using enable_if_neq_comparable = decltype(concepts::type_require<T>(std::declval<I>() != std::declval<S>()));

template<typename I, typename S>
struct iterator_traits<I,
                       S,
                       enable_if_t<std::is_base_of<iterator_facade<I>, I>::value,
                                   enable_if_neq_comparable<I, S, typename I::iterator_category>>>
{};

/// A non owning range CRTP base
template<typename Derived, typename I, typename S = I> struct view_facade
{

  private:
    using it_traits = iterator_traits<I, S>;

    template<typename D = Derived, typename = enable_if_t<std::is_same<D, Derived>::value>>
    auto derived() -> decltype(*static_cast<D*>(this))
    {
        return *static_cast<Derived*>(this);
    }
    template<typename D = Derived, typename = enable_if_t<std::is_same<D, Derived>::value>>
    auto derived() const -> decltype(*static_cast<const D*>(this))
    {
        return *static_cast<const Derived*>(this);
    }

  public:
    using iterator        = I;
    using const_iterator  = iterator;
    using sentinel        = S;
    using reference       = typename it_traits::reference;
    using const_reference = typename it_traits::reference;
    using element_type    = typename it_traits::value_type;
    using difference_type = typename it_traits::difference_type;
    using size_type       = typename it_traits::size_type;

    bool empty() const
    {
        const auto first = derived().begin();
        return !bound_check(first, first, derived().end());
    }

    template<typename D = Derived>
    auto size() const noexcept -> decltype(it_traits::size(derived<D>().begin(), derived<D>().end()))
    {
        return it_traits::size(derived().begin(), derived().end());
    }

    iterator cbegin() const noexcept { return derived().begin(); }
    sentinel cend() const noexcept { return derived().end(); }

    template<typename D = Derived>
    auto rbegin() const
      -> enable_if_t<iterator_traits<S, I>::is_bidirectional, std::reverse_iterator<decltype(derived<D>().end())>>
    {
        return std::reverse_iterator<S>{derived().end()};
    }
    template<typename D = Derived>
    auto rend() const
      -> enable_if_t<iterator_traits<S, I>::is_bidirectional, std::reverse_iterator<decltype(derived<D>().begin())>>
    {
        return std::reverse_iterator<I>{derived().begin()};
    }

    template<typename D = Derived> auto crbegin() const noexcept -> decltype(derived<D>().rbegin())
    {
        return derived().rbegin();
    }
    template<typename D = Derived> auto crend() const noexcept -> decltype(derived<D>().rend())
    {
        return derived().rend();
    }

    template<typename D = Derived>
    auto front() const noexcept -> decltype(it_traits::front(derived<D>().begin(), derived<D>().end()))
    {
        return it_traits::front(derived().begin(), derived().end());
    }
    template<typename D = Derived>
    auto back() const noexcept -> decltype(it_traits::back(derived<D>().begin(), derived<D>().end()))
    {
        return it_traits::front(derived().begin(), derived().end());
    }

    template<typename D = Derived>
    auto operator[](size_type idx) const -> decltype(it_traits::at(idx, derived<D>().begin(), derived<D>().end()))
    {
        return it_traits::at(idx, derived().begin(), derived().end());
    }

    template<typename D = Derived, typename = enable_if_t<std::is_same<D, Derived>::value>>
    friend std::ostream& operator<<(std::ostream& sout, const Derived& rng)
    { // Courtesy of range-v3
        sout << '[';
        auto       it = rng.begin();
        auto const e  = rng.end();
        if (it != e) {
            for (;;) {
                sout << *it;
                if (++it == e) break;
                sout << ',';
            }
        }
        sout << ']';
        return sout;
    }
};

template<typename R, typename Traits>
auto
operator<<(std::basic_ostream<char, Traits>& os, const R& s)
  -> decltype(os.write(static_cast<const char*>(begin(s)), size(s)))
{
    return os.write(begin(s), size(s));
}

template<typename R,
         typename Char   = char,
         typename Traits = std::char_traits<Char>,
         typename Alloc  = std::allocator<Char>>
auto
to_string(const R& r) -> decltype(std::basic_string<Char, Traits, Alloc>(begin(r), size(r)))
{
    return std::basic_string<Char, Traits, Alloc>(begin(r), size(r));
}

/// Find specialization for char like types
template<typename T, typename = enable_if_t<sizeof(T) == 1>>
inline T*
find(T* first, const T* last, T v)
{
    assume(first <= last, "swapped iterators");
    auto p = static_cast<T*>(memchr(first, reinterpret_cast<const char&>(v), size_t(last - first)));
    return likely(p != nullptr) ? p : const_cast<T*>(last);
}

template<typename T, typename = enable_if_t<sizeof(T) == 1>>
inline const T*
find(const T* first, const T* last, T v)
{
    assume(first <= last, "swapped iterators");
    const T* p = static_cast<const T*>(memchr(first, reinterpret_cast<const char&>(v), size_t(last - first)));
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

#endif // GATBL_RANGES_HPP
