#ifndef REVERSE_RANGE_HPP
#define REVERSE_RANGE_HPP

#include "gatbl/utils/ranges.hpp"

namespace gatbl {

template<typename Repr> struct reverse_range
{
    using element_type           = value_t<Repr>;
    using reverse_iterator       = iterator_t<Repr>;
    using const_reverse_iterator = iterator_t<const Repr>;
    using iterator               = std::reverse_iterator<reverse_iterator>;
    using const_iterator         = std::reverse_iterator<const_reverse_iterator>;

    // public for brace initialization
    Repr _data;

    constexpr auto size() const noexcept -> decltype(size(this->_data)) { return size(this->_data); }
    constexpr auto empty() const noexcept -> decltype(empty(this->_data)) { return empty(this->_data); }

    friend constexpr const_iterator begin(const reverse_range& r) noexcept { return const_iterator(end(r._data)); }
    friend constexpr const_iterator end(const reverse_range& r) noexcept { return const_iterator(begin(r._data)); }
    friend constexpr iterator       begin(reverse_range& r) noexcept { return iterator(end(r._data)); }
    friend constexpr iterator       end(reverse_range& r) noexcept { return iterator(begin(r._data)); }

    friend constexpr const_reverse_iterator rbegin(const reverse_range& r) noexcept { return begin(r._data); }
    friend constexpr const_reverse_iterator rend(const reverse_range& r) noexcept { return end(r._data); }
    friend constexpr reverse_iterator       rbegin(reverse_range& r) noexcept { return begin(r._data); }
    friend constexpr reverse_iterator       rend(reverse_range& r) noexcept { return end(r._data); }

    constexpr element_type&       front() noexcept { return *begin(*this); }
    constexpr const element_type& front() const noexcept { return *begin(*this); }
    constexpr element_type&       back() noexcept
    {
        assume(!empty());
        return *(end(*this) - 1);
    }
    constexpr const element_type& back() const noexcept
    {
        assume(!empty());
        return *(end(*this) - 1);
    }

    const element_type& operator[](size_t i) const noexcept
    {
        assume(i < size());
        return *(begin(*this) + i);
    }

    element_type& operator[](size_t i) noexcept
    {
        assume(i < size());
        return *(begin(*this) + i);
    }
};

}

#endif // REVERSE_RANGE_HPP
