#ifndef CONCEPTS_H
#define CONCEPTS_H

#include <utility>

namespace gatbl {
/** While waiting for real concepts adoption in compiler this library proposeto use template function signature to
 * express requirement on type.
 *
 * For example write a function that check for members :
 * template<typename T> auto plop(T& x) -> decltype(type_require<return_t>(x.plip())) { ... }
 * or
 * template<typename T> auto plop(T& x) -> decltype(value_require(c.compute_return(), x.plip())) { ... }
 *
 * Different headers in gatbl extends the concepts namespace (e.g see ranges.hpp)
 */
namespace concepts {

///
template<typename T, typename... X>
auto
type_require(X&&...) -> T;

template<typename X, typename... Xs>
constexpr inline auto
value_require(X&& x, Xs&&...) -> decltype(std::forward<X>(x))
{
    return std::forward<X>(x);
}

template<typename To, typename From>
auto
ImplicitConvertible(From&& x) -> decltype(type_require<To, To>(std::forward<From>(x)));
// template<typename To, typename From> auto Convertible(From&& x) -> decltype(convert<To>(std::forward<From>(x)));

template<typename To, typename From = To>
auto
Assignable(To& to, From&& from) -> decltype(to = std::forward<From>(from));

/// Tag types from any type (for function argument tags mainly, named after Proxy in Haskell)
/// eg: template<typename R> R func(..., proxy_t<R>={}); allows to chose overload's return type with func(...,
/// proxy_t<return_t>{})

} // namespace concepts

template<typename... T> struct proxy_t
{};

} // namespace gatbl

#endif // CONCEPTS_H
