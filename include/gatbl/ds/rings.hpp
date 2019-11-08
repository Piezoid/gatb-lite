#ifndef RINGS_HPP
#define RINGS_HPP

#include <memory>
#include <limits>

#include "gatbl/common.hpp"
#include "gatbl/sys/memory.hpp"

namespace gatbl {

template<typename T> struct wrapper
{
    wrapper() = default;
    wrapper(T v)
      : val(v)
    {}
    wrapper(const wrapper&) = default;
    wrapper& operator=(const wrapper&) = default;
    wrapper& operator                  =(const T& v)
    {
        val = v;
        return *this;
    }
      operator const T&() const { return val; }
      operator T&() { return val; }
    T val;
};

template<typename elem_t, size_t log2size = 8, typename idx_t = uint_fast8_t, bool heap = true> struct ring_deque
{
    static constexpr size_t capacity = 1ull << log2size;
    static_assert(capacity - 1 <= std::numeric_limits<idx_t>::max(), "idx_t too short");

    ring_deque()                  = default;
    ring_deque(const ring_deque&) = delete;
    ring_deque(ring_deque&&)      = default;

    ~ring_deque()
    {
        const auto end = _e;
        for (idx_t i = _b; i != end; ++i)
            _arr.pop_at(i);
    }

    template<typename... Args> elem_t& emplace_front(Args&&... args)
    {
        check_not_full();
        return _arr.emplace_at(mask & --_b, std::forward<Args>(args)...);
    }
    template<typename... Args> elem_t& emplace_back(Args&&... args)
    {
        check_not_full();
        return _arr.emplace_at(mask & _e++, std::forward<Args>(args)...);
    }
    elem_t& front()
    {
        check_not_empty();
        return _arr[mask & _b];
    }
    elem_t& back()
    {
        check_not_empty();
        return _arr[mask & idx_t(_e - 1)];
    }
    const elem_t& front() const
    {
        check_not_empty();
        return _arr[mask & _b];
    }
    const elem_t& back() const
    {
        check_not_empty();
        return _arr[mask & idx_t(_e - 1)];
    }
    elem_t pop_front()
    {
        check_not_empty();
        return _arr.pop_at(mask & _b++);
    }
    elem_t pop_back()
    {
        check_not_empty();
        return _arr.pop_at(mask & --_e);
    }
    void  clear() { _e = _b; }
    bool  empty() const { return _b == _e; }
    bool  full() const { return size() == capacity; }
    idx_t size() const
    {
        assume(idx_t(_e - _b) <= capacity, "Unwrapped indices");
        return _e - _b;
    }

  private:
    void check_not_full() const { assume(!full(), "Ring full"); }
    void check_not_empty() const { assume(!empty(), "Ring empty"); }

    static constexpr idx_t                              mask = capacity - 1;
    idx_t                                               _b   = 0;
    idx_t                                               _e   = 0;
    unsafe::uninitialized_array<elem_t, capacity, heap> _arr;
};

} /* namespace gatbl */

#endif // RINGS_HPP
