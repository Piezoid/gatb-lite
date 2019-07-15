#ifndef NO_CONVERSION_HPP
#define NO_CONVERSION_HPP

#include <type_traits>
#include <utility>

namespace gatbl {

template<typename T, typename U, typename Enabled = void> struct is_same_cvref : std::false_type
{};

template<typename T, typename E> struct is_same_cvref<T, T, E> : std::true_type
{
    using enable_t = E;
};

template<typename T, typename E> struct is_same_cvref<T, T&, E> : std::true_type
{
    using enable_t = E;
};

template<typename T, typename E> struct is_same_cvref<T, T&&, E> : std::true_type
{
    using enable_t = E;
};

template<typename T, typename E> struct is_same_cvref<T, const T, E> : std::true_type
{
    using enable_t = E;
};

template<typename T, typename E> struct is_same_cvref<T, const T&, E> : std::true_type
{
    using enable_t = E;
};

template<typename T, typename E> struct is_same_cvref<T, const T&&, E> : std::true_type
{
    using enable_t = E;
};

template<typename T, typename U> constexpr bool is_same_cvref_v = is_same_cvref<T, U>::value;

template<typename T, typename U, typename Enabled = void>
using if_is_same_cvref = typename is_same_cvref<T, U, Enabled>::enable_t;

/// A wrapper for arguments that disable implicit conversions
/// For example `void foo(no_conversion<uint16_t>)` can't be invoked as foo(uint8_t(0)).
template<typename T> struct no_conversion
{
    template<typename U, typename = if_is_same_cvref<T, U>>
    no_conversion(U&& v) // Implicit constructor with wildcard type argument for catching all conversion attempts
      : v_(std::forward<U>(v))
    {}
    operator T() const { return v_; }

  private:
    T v_;
};

} // namespace gatbl

#endif // NO_CONVERSION_HPP
