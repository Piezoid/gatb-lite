#ifndef UNIQUE_RANGE_HPP
#define UNIQUE_RANGE_HPP

#include <memory>

namespace gatbl {

template<typename T, typename D = std::default_delete<T[]>> struct unique_range : private D
{

    using value_type      = T;
    using reference       = value_type&;
    using const_reference = const value_type&;
    using iterator        = value_type*;
    using const_iterator  = const value_type*;
    using sentinel        = const_iterator;

    constexpr unique_range() = default;

    template<typename _D = D>
    unique_range(T* ptr, size_t size, _D&& d = D{}) noexcept

      : D(std::forward<_D>(d))
      , _ptr(ptr)
      , _size(size)
    {}
    ~unique_range() { reset(); }
    unique_range(const unique_range&) = delete;
    unique_range(unique_range&& from) noexcept
    {
        std::swap(from._ptr, _ptr);
        std::swap(from._size, _size);
    }
    unique_range& operator=(const unique_range&) = delete;
    unique_range& operator                       =(unique_range&& from) noexcept
    {
        std::swap(from._ptr, _ptr);
        std::swap(from._size, _size);
        return *this;
    }

    size_t size() const noexcept { return _ptr ? _size : 0; }
    bool   empty() const noexcept { return _ptr && _size; }
           operator bool() const noexcept { return empty(); }

    const_iterator begin() const noexcept { return _ptr; }
    const_iterator end() const noexcept { return begin() + size(); }

    operator T*() const { return begin(); }

    const_reference operator[](size_t i) const noexcept
    {
        assume(i < size(), "out of bound i=%lu < size()=%lu", i, size());
        return *(begin() + i);
    }

    void reset(std::nullptr_t = nullptr)
    {
        T*     ptr  = nullptr;
        size_t size = 0;

        std::swap(_ptr, ptr);
        std::swap(_size, size);

        if (ptr != nullptr) call_deleter(*static_cast<D*>(this), ptr, size);
    }

  protected:
    template<typename _D> static auto call_deleter(_D& d, T* ptr, size_t size) -> decltype(d(ptr, size))
    {
        return d(ptr, size);
    }
    template<typename _D> static auto call_deleter(_D& d, T* ptr, size_t) -> decltype(d(ptr)) { return d(ptr); }

    T*     _ptr  = nullptr;
    size_t _size = 0;
};

} // namespace gatbl

#endif // UNIQUE_RANGE_HPP
