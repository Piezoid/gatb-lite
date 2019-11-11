#ifndef GATBL_NUCLEOTIDE_HPP
#define GATBL_NUCLEOTIDE_HPP

#include "gatbl/common.hpp"
#include "gatbl/utils/compatibility.hpp"

namespace gatbl {

namespace details {
template<typename T, typename = typename std::is_enum<T>::type> struct _to_underlying
{
    using type = T;
    constexpr static T call(T v) noexcept { return v; }
};

template<typename E> struct _to_underlying<E, std::true_type>
{
    using type = typename std::underlying_type<E>::type;
    constexpr static type call(E v) noexcept { return static_cast<type>(v); }
};

/// For the sake of bit_vector, bool is considered as an enum
template<> struct _to_underlying<bool, std::false_type>
{
    using type = uint8_t;
    constexpr static type call(bool v) noexcept { return static_cast<type>(v); }
};

}

template<typename T> using underlying_type = typename details::_to_underlying<T>::type;

/// If `T` is an enum or bool, returns the underlying type, such that T(underlying_type(x)) is the identity function.
/// Otherwise this is the identity function.
template<typename T>
static constexpr underlying_type<T>
to_underlying(T v) noexcept
{
    return details::_to_underlying<T>::call(v);
}

using nucint_t = uint_fast8_t;
enum class nuc_t : nucint_t { A = 0, C, T, G, N = 255 };

/// Same as to_underlying(nuc_t x) but performs domain check in debug builds
template<typename T = nucint_t>
static inline enable_if_t<T(0) < T(-1), T>
as_integer(nuc_t n)
{
    auto i = static_cast<nucint_t>(n);
    assume(i < 4, "invalid nucleotide %u", unsigned(i));
    return i;
}

/// Same as nuc_t(nucint_t x) but performs domain check in debug builds
template<typename T>
static inline enable_if_t<T(0) < T(-1), nuc_t>
as_nuc(T n)
{
    assume(n < 4 || n == as_integer<T>(nuc_t::N), "invalid nucleotide %u", unsigned(n));
    return nuc_t(n);
}

inline nucint_t
complement(nucint_t nuc)
{
    return nuc ^ nucint_t(2u);
}

inline nuc_t
complement(nuc_t nuc)
{
    return as_nuc(complement(as_integer(nuc)));
}

inline std::ostream&
operator<<(std::ostream& o, nuc_t n)
{
    return o.put("ACTG"[as_integer(n)]);
}

}

#endif // GATBL_NUCLEOTIDE_HPP
