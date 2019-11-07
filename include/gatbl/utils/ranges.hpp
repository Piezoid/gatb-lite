#ifndef GATBL_RANGES_HPP
#define GATBL_RANGES_HPP

#include <cstring>
#include "gatbl/common.hpp"
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

namespace concepts {

template<typename R>
auto
Range(R&& r) -> decltype(type_require<R>(begin(r) != end(r), ++begin(r)));

template<typename R>
auto
FiniteRange(R&& r) -> decltype(type_require<R>(size(r), begin(r) != end(r), ++begin(r)));

template<typename T, typename R>
auto
RangeOf(R&& r) -> decltype(value_require<T&>(*++begin(r), begin(r) != end(r)));

template<typename T, typename R>
auto
FiniteRangeOf(R&& r) -> decltype(value_require<T&>(*++begin(r), begin(r) != end(r), size(r)));

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

template<typename L, typename R>
constexpr static auto
compare(const L& lhs, const R& rhs) noexcept -> decltype(lhs - rhs)
{
    return lhs - rhs;
}

/// Comparable CRTP base, user must define either `signed operator-(L, R)` or `signed compare(L, R)`.
/// Otherwise all of the following must be defined:
/// - operator<(L, R)
/// - operator<(R, L)
/// - operator==(R, L)
template<typename L, typename R = L> struct comparable_facade
{
  private:
    const L& derived() const noexcept { return static_cast<const L&>(*this); }

  public:
    bool operator==(const R& other) const { return compare(derived(), other) == 0; }
    bool operator<(const R& other) const { return compare(derived(), other) < 0; }
    bool operator>(const R& other) const { return compare(derived(), other) > 0; }
    bool operator!=(const R& other) const { return !(derived() == other); }
    bool operator>=(const R& other) const { return !(derived() < other); }
    bool operator<=(const R& other) const { return !(derived() > other); }
};

/// Iterator CRTP base, similar to std::iterator, but with convenience methods
template<typename D,
         typename Category,
         typename ValueTy,
         typename Reference = ValueTy&,
         typename Pointer   = ValueTy*,
         typename Distance  = ptrdiff_t,
         typename Sentinel  = D>
struct iterator_facade;

template<typename D, typename ValueTy, typename Reference, typename Pointer, typename Distance, typename Sentinel>
struct iterator_facade<D, gatbl::contiguous_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>

  : public iterator_facade<D, std::random_access_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
{
    using iterator_category = gatbl::contiguous_iterator_tag;
};

template<typename D, typename ValueTy, typename Reference, typename Pointer, typename Distance, typename Sentinel>
struct iterator_facade<D, std::random_access_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
  : iterator_facade<D, std::bidirectional_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
  , comparable_facade<D, Sentinel>
{
  protected:
    using iterator_facade<D, std::bidirectional_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>::derived;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using comparable_facade<D, Sentinel>::operator==;
    using comparable_facade<D, Sentinel>::operator!=;

    CPP14_CONSTEXPR D& operator+=(Distance x) { return derived().operator-=(-x); }

    CPP14_CONSTEXPR D& operator-=(Distance x) { return derived().operator+=(-x); }

    CPP14_CONSTEXPR D operator+(Distance x) const
    {
        D   tmp = derived();
        tmp.operator+=(x);
        return tmp;
    }

    friend constexpr D operator+(Distance x, D const& it) { return it + x; }

    CPP14_CONSTEXPR D operator-(Distance x) const
    {
        D   tmp = derived();
        tmp.operator-=(x);
        return tmp;
    }

    constexpr Reference operator[](Distance x) const { return *(derived() + x); }
};

template<typename D, typename ValueTy, typename Reference, typename Pointer, typename Distance, typename Sentinel>
struct iterator_facade<D, std::bidirectional_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
  : iterator_facade<D, std::forward_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
{
    using iterator_category = std::bidirectional_iterator_tag;

    D& operator--() { return this->derived().operator-=(1); }

    CPP14_CONSTEXPR D operator--(int)
    {
        auto& ref = static_cast<D&>(*this);
        D     tmp = ref;
        ref.  operator--();
        return tmp;
    }
};

template<typename D, typename ValueTy, typename Reference, typename Pointer, typename Distance, typename Sentinel>
struct iterator_facade<D, std::forward_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
  : iterator_facade<D, std::input_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
{
  public:
    using iterator_category = std::forward_iterator_tag;
    using iterator_facade<D, std::input_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>::operator++;

    CPP14_CONSTEXPR D operator++(int)
    {
        auto& ref = static_cast<D&>(*this);
        D     tmp = ref;
        ref.  operator++();
        return tmp;
    }
};

template<typename D, typename ValueTy, typename Reference, typename Pointer, typename Distance, typename Sentinel>
struct iterator_facade<D, std::input_iterator_tag, ValueTy, Reference, Pointer, Distance, Sentinel>
{
  protected:
    D&       derived() noexcept { return static_cast<D&>(*this); }
    const D& derived() const noexcept { return static_cast<const D&>(*this); }

  public:
    using iterator_category = std::input_iterator_tag;
    using difference_type   = Distance;
    using pointer           = Pointer;
    using reference         = Reference;
    using value_type        = ValueTy;

    constexpr bool operator!=(const D& other) const { return !(derived() == other); }

    constexpr bool operator==(const D& other) const { return !(derived() != other); }

    CPP14_CONSTEXPR D& operator++() { return derived().operator+=(1u); }
};

/// A non owning range CRTP base
template<typename Derived,
         typename I,
         typename CI    = I,
         typename S     = I,
         typename CS    = CI,
         typename ItCat = typename std::iterator_traits<I>::iterator_category>
struct view_facade;

template<typename Derived, typename I, typename CI, typename S, typename CS>
struct view_facade<Derived, I, CI, S, CS, gatbl::contiguous_iterator_tag>
  : view_facade<Derived, I, CI, S, CS, std::random_access_iterator_tag>
{
  private:
    using base = view_facade<Derived, I, CI, S, CS, std::input_iterator_tag>;

  public:
    using typename base::size_type;
    using typename base::value_type;
    using pointer       = value_type*;
    using const_pointer = const value_type*;

    const pointer* data() const { return std::addressof(*this->derived().begin()); }
    const_pointer  data() { return std::addressof(*this->derived().begin()); }
};

template<typename Derived, typename I, typename CI, typename S, typename CS>
struct view_facade<Derived, I, CI, S, CS, std::random_access_iterator_tag>
  : view_facade<Derived, I, CI, S, CS, std::bidirectional_iterator_tag>
{
  private:
    using base = view_facade<Derived, I, CI, S, CS, std::bidirectional_iterator_tag>;

  public:
    using typename base::const_reference;
    using typename base::reference;
    using typename base::size_type;

    size_type size() const { return size_type(this->derived().end() - this->derived().begin()); }

    const_reference at(size_type i) const
    {
        assume(i < this->derived().size(), "out of bound indice");
        return *(this->derived().begin() + i);
    }

    reference at(size_type i)
    {
        assume(i < this->derived().size(), "out of bound indice");
        return *(this->derived().begin() + i);
    }
};

template<typename Derived, typename I, typename CI, typename S, typename CS>
struct view_facade<Derived, I, CI, S, CS, std::bidirectional_iterator_tag>
  : view_facade<Derived, I, CI, S, CS, std::forward_iterator_tag>
{
  private:
    using base = view_facade<Derived, I, CI, S, CS, std::forward_iterator_tag>;

  public:
    using typename base::const_reference;
    using typename base::reference;
    using reverse_iterator       = std::reverse_iterator<I>;
    using const_reverse_iterator = std::reverse_iterator<I>;

    reverse_iterator       rbegin() { return std::reverse_iterator<S>{this->derived().end()}; }
    reverse_iterator       rend() { return std::reverse_iterator<I>{this->derived().begin()}; }
    const_reverse_iterator rbegin() const { return std::reverse_iterator<S>{this->derived().end()}; }
    const_reverse_iterator rend() const { return std::reverse_iterator<I>{this->derived().begin()}; }
    const_reverse_iterator crbegin() const noexcept { return this->derived().rbegin(); }
    const_reverse_iterator crend() const noexcept { return this->derived().rend(); }

    const_reference back() const
    {
        this->check_not_empty();
        auto it = this->derived().end();
        --it;
        return *it;
    }

    reference back()
    {
        this->check_not_empty();
        auto it = this->derived().end();
        --it;
        return *it;
    }
};

template<typename Derived, typename I, typename CI, typename S, typename CS>
struct view_facade<Derived, I, CI, S, CS, std::forward_iterator_tag>
  : view_facade<Derived, I, CI, S, CS, std::input_iterator_tag>
{
  private:
    using base = view_facade<Derived, I, CI, S, CS, std::input_iterator_tag>;

  public:
    using typename base::const_reference;
    using typename base::reference;
    using typename base::size_type;

    I cbegin() const noexcept { return this->derived().begin(); }
    S cend() const noexcept { return this->derived().end(); }

    size_type size() const { return size_type(std::distance(this->derived().begin(), this->derivend().end())); }

    const_reference at(size_type i) const
    {
        using std::advance;
        auto it = this->derived().begin();
        advance(it, i);
        assert(it != this->derived().begin(), "out of bound indice");
        return *it;
    }

    reference at(size_type i)
    {
        using std::advance;
        auto it = this->derived().begin();
        advance(it, i);
        assert(it != this->derived().begin(), "out of bound indice");
        return *it;
    }

    const_reference operator[](size_type i) const { return this->derived().at(i); }
    reference       operator[](size_type i) { return this->derived().at(i); }

    const_reference back() const
    {
        this->check_not_empty();
        return this->derived().at(this->derived().size() - 1u);
    }

    reference back()
    {
        this->check_not_empty();
        return this->derived().at(this->derived().size() - 1u);
    }
};

template<typename Derived, typename I, typename CI, typename S, typename CS>
struct view_facade<Derived, I, CI, S, CS, std::input_iterator_tag>
{
  protected:
    using it_traits = std::iterator_traits<I>;

    const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    Derived&       derived() noexcept { return static_cast<Derived&>(*this); }

    void check_not_empty() const { assert(not empty(), "empty range"); }

  public:
    using iterator        = I;
    using const_iterator  = CI;
    using sentinel        = S;
    using const_sentinel  = CS;
    using reference       = typename it_traits::reference;
    using const_reference = typename std::iterator_traits<CI>::reference;
    using value_type      = typename it_traits::value_type;
    using difference_type = typename it_traits::difference_type;
    using size_type       = make_unsigned_t<difference_type>;

    bool empty() const noexcept { return derived().begin() == derived().end(); }

    template<typename D = Derived> const_reference front() const noexcept
    {
        check_not_empty();
        return derived().begin();
    }

    template<typename D = Derived> reference front() noexcept
    {
        check_not_empty();
        return derived().begin();
    }

    const_iterator cbegin() const noexcept { return derived().begin(); }
    const_sentinel cend() const noexcept { return derived().end(); }

    //    friend auto operator<<(std::ostream& sout, const Derived& rng) -> decltype(sout <<
    //    std::declval<value_type>()) { // Courtesy of range-v3
    //        sout << '[';
    //        auto       it = rng.derived().begin();
    //        auto const e  = rng.derived().end();
    //        if (it != e) {
    //            for (;;) {
    //                sout << *it;
    //                if (++it == e) break;
    //                sout << ',';
    //            }
    //        }
    //        sout << ']';
    //        return sout;
    //    }
};

// String conversion
template<typename R, typename CharT = char, typename Traits = std::char_traits<CharT>>
auto
operator<<(std::basic_ostream<CharT, Traits>& os, const R& s)
  -> decltype(os.write(static_cast<const CharT*>(begin(s)), size(s)))
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

using std::find;

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

}

#endif // GATBL_RANGES_HPP
