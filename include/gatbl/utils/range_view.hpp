#ifndef RANGE_VIEW_H
#define RANGE_VIEW_H

#include "gatbl/utils/empty_base.hpp"
#include "gatbl/utils/ranges.hpp"
#include "gatbl/utils/concepts.h"

namespace gatbl {

/// Construct a range wrapper giver an isomorphism functor (with one overload of operator() for each side of the
/// bijection)
template<typename Rep, typename Iso> class range_view : public utils::empty_base<Rep, Iso>
{
    using base      = utils::empty_base<Rep, Iso>;
    using _iterator = iterator_t<Rep>;

    template<typename T> class wreference : protected utils::empty_base<T, Iso>
    {
        using base = utils::empty_base<T, Iso>;

      protected:
        using base::base;

      public:
        operator auto() const { return this->tag()(this->value()); }

        template<typename V>
        auto operator=(V&& v)
          -> decltype(concepts::type_require<wreference&>(this->value() = this->tag()(std::forward<V>(v))))
        {
            this->value() = this->tag()(std::forward<V>(v));
            return *this;
        }
    };

    template<typename T> class witerator : protected utils::empty_base<T, Iso>
    {
        using _value_type = decltype(*std::declval<T>());
        using base        = utils::empty_base<T, Iso>;
        using base::base;

      public:
        witerator& operator++()
        {
            ++this->value();
            return *this;
        }
        bool                    operator!=(witerator& other) { return this->value() != other.value(); }
        wreference<_value_type> operator*() { return wreference<_value_type>(*(this->value()), this->tag()); }
    };

  public:
    using base::base;

    using iterator               = witerator<iterator_t<Rep>>;
    using const_iterator         = witerator<iterator_t<const Rep>>;
    using sentinel               = witerator<sentinel_t<Rep>>;
    using const_sentinel         = witerator<sentinel_t<const Rep>>;
    using reverse_iterator       = witerator<reverse_iterator_t<Rep>>;
    using reverse_const_iterator = witerator<reverse_iterator_t<const Rep>>;
    using reverse_sentinel       = witerator<reverse_sentinel_t<Rep>>;
    using reverse_const_sentinel = witerator<reverse_sentinel_t<const Rep>>;

    friend constexpr iterator               begin(range_view& rv) { return {begin(rv.value()), rv.tag()}; }
    friend constexpr const_iterator         begin(const range_view& rv) { return {begin(rv.value()), rv.tag()}; }
    friend constexpr sentinel               end(range_view& rv) { return {end(rv.value()), rv.tag()}; }
    friend constexpr const_sentinel         end(const range_view& rv) { return {end(rv.value()), rv.tag()}; }
    friend constexpr reverse_iterator       rbegin(range_view& rv) { return {rbegin(rv.value()), rv.tag()}; }
    friend constexpr reverse_const_iterator rbegin(const range_view& rv) { return {rbegin(rv.value()), rv.tag()}; }
    friend constexpr reverse_sentinel       rend(range_view& rv) { return {rend(rv.value()), rv.tag()}; }
    friend constexpr reverse_const_sentinel rend(const range_view& rv) { return {rend(rv.value()), rv.tag()}; }
};

}

#endif // RANGE_VIEW_H
