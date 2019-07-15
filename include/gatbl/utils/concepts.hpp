#ifndef CONCEPTS_H
#define CONCEPTS_H

#include <cstddef>
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
T
type_require(X&&...);

template<typename X, typename... Xs>
constexpr inline auto
value_require(X&& x, Xs&&...) -> decltype(std::forward<X>(x));

template<typename To, typename From>
auto
ImplicitConvertible(From&& x) -> decltype(type_require<To, To>(std::forward<From>(x)));

template<typename To, typename From = To>
auto
Assignable(To& to, From&& from) -> decltype(to = std::forward<From>(from));

/** Allow to match on array types for SFINAE
 * eg: std::unique_ptr<typename extent_kind<X>::unknown_bound_array[]>
 * is valid for int[], but not int or int[4]
 */
template<typename T> struct extent_kind
{
    using base_type     = T;
    using single_object = T;
};

template<class T> struct extent_kind<T[]>
{
    using base_type           = T;
    using unknown_bound_array = T;
    using array               = T;
};

template<class T, size_t N> struct extent_kind<T[N]>
{
    using base_type         = T;
    using known_bound_array = T;
    using array             = T;
};

/// Tag types from any type (for function argument tags mainly, named after Proxy in Haskell)
/// eg: template<typename R> R func(..., proxy_t<R>={}); allows to chose overload's return type with func(...,
/// proxy_t<return_t>{})
template<typename... T> struct proxy_t
{};

} // namespace concepts
} // namespace gatbl

#endif // CONCEPTS_H
