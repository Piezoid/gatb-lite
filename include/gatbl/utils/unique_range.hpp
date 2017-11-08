#ifndef UNIQUE_RANGE_HPP
#define UNIQUE_RANGE_HPP

#include <memory>
#include "gatbl/utils/empty_base.hpp"

namespace gatbl {

template<typename T, typename D = std::default_delete<T[]>> struct unique_range
{

    using element_type    = T;
    using reference       = element_type&;
    using const_reference = const element_type&;
    using iterator        = element_type*;
    using const_iterator  = const element_type*;
    using sentinel        = const_iterator;

    constexpr unique_range() = default;

    template<typename _D = D>
    unique_range(T* ptr, size_t size, _D&& d = D{}) noexcept
      : _ptr(ptr, std::forward<D>(d))
      , _size(size)
    {}
    unique_range(const unique_range&) = delete;
    unique_range(unique_range&& from) noexcept
      : unique_range()
    {
        std::swap(from._ptr, _ptr);
        std::swap(from._size, _size);
    }
    unique_range& operator=(const unique_range&) = delete;
    unique_range& operator                       =(unique_range&& from) noexcept
    {
        reset();
        std::swap(from._ptr, _ptr);
        std::swap(from._size, _size);
        return *this;
    }

    size_t size() const noexcept { return _ptr.value() ? _size : 0; }
    bool   empty() const noexcept { return _ptr.value() && _size; }

    iterator       begin() noexcept { return _ptr.value(); }
    const_iterator begin() const noexcept { return _ptr.value(); }
    iterator       end() noexcept { return begin() + size(); }
    const_iterator end() const noexcept { return begin() + size(); }

    reference operator[](size_t i) noexcept
    {
        assume(i < size());
        return *(begin() + i);
    }

    const_reference operator[](size_t i) const noexcept
    {
        assume(i < size());
        return *(begin() + i);
    }

    void reset(std::nullptr_t = nullptr)
    {
        if (!empty()) { destruct(_ptr.tag()); }
        _ptr  = nullptr;
        _size = 0;
    }

  protected:
    // FIXME: would like to enable this for default_delete<T[]>
    //    auto destruct(D& d) -> decltype(d(this->begin()))
    //    {
    //        utils::empty_base<T*, D> ptr = { nullptr, {} };
    //        std::swap(_ptr, ptr);
    //        _size = 0;
    //        return ptr.tag()(ptr.value());
    //    }

    void destruct(D& d /*, std::nullptr_t = nullptr*/) //-> decltype(d(this->begin(), this->size()))
    {
        utils::empty_base<T*, D> ptr = {nullptr, {}};
        std::swap(_ptr, ptr);
        size_t size = 0;
        std::swap(_size, size);
        return ptr.tag()(ptr.value(), size);
    }

    utils::empty_base<T*, D> _ptr  = {nullptr, {}};
    size_t                   _size = 0;
};

} // namespace gatbl

#endif // UNIQUE_RANGE_HPP
