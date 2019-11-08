#ifndef MEMORY_HPP
#define MEMORY_HPP

#include <memory>

#include "gatbl/common.hpp"
#include "gatbl/utils/compatibility.hpp"

// Purpose:
// Before C++17's new doesn't correctly align over-aligned types, so we use posix_memalign and a free() based
// deallocator

namespace gatbl {
namespace details {

template<typename T> struct free_delete
{
    void operator()(T* ptr) const
    {
        assume(ptr != nullptr, "free_delete(nullptr)");
        ptr->~T();
        free(ptr);
    }
};

template<typename T> struct free_delete<T[]>
{
    // Memory allocated for arrays of non-trivially-destructible objects have a header before the first object
    static constexpr size_t header_size
      = std::is_trivially_destructible<T>::value ? 0 : (sizeof(T) > sizeof(size_t) ? sizeof(T) : sizeof(size_t));

    /// Allows instances to be used as a deleter
    void operator()(T* ptr) const
    {
        assume(ptr != nullptr, "free_delete(nullptr)");
        CPP17_IF_CONSTEXPR(header_size > 0)
        { // non-trivially-destructible objects
            // Retreive the size of the array stored before the first element
            for (T* it = ptr + *(reinterpret_cast<const size_t*>(ptr) - 1); it-- > ptr;)
                it->~T();
        }

        // Retreive the original pointer from the malloc/aligned_alloc allocation
        free(reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(ptr) - header_size));
    }

    static T* allocate(size_t n)
    {
        const size_t alloc_bytes = n * sizeof(T) + header_size;
        void*        void_ptr;

        constexpr size_t align = alignof(T) > sizeof(void*) ? alignof(T) : sizeof(void*);
        if (posix_memalign(&void_ptr, align, alloc_bytes) != 0) { throw std::bad_alloc(); }
        auto first_object = reinterpret_cast<uintptr_t>(void_ptr) + header_size;

        // Store the array size before the first object
        *(reinterpret_cast<size_t*>(first_object) - 1) = n;

        return reinterpret_cast<T*>(first_object);
    }
};

} /* namspace details */

/** A unique pointer to a malloc/aligned_alloced object */
template<typename T> using unique_malloc_ptr = std::unique_ptr<T, details::free_delete<T>>;

namespace unsafe {

/** Allocates memory for an object of type T, respecting overalignement requirement
 * The object are not initialized !
 * Because the destructor will eventually be called, the user must initialize the object if it is not trivially
 * destructible (eg. containing pointers of any sort).
 */
template<typename T>
unique_malloc_ptr<typename concepts::extent_kind<T>::single_object>
make_unique_aligned_uninitialized()
{
    void*            void_ptr;
    constexpr size_t align = alignof(T) > sizeof(void*) ? alignof(T) : sizeof(void*);
    if (posix_memalign(&void_ptr, align, sizeof(T)) != 0) { throw std::bad_alloc(); }
    return unique_malloc_ptr<T>(reinterpret_cast<T*>(void_ptr));
}

/** Allocates a slice of memory for n objects of type T, respecting overalignement requirement
 * The objects are not initialized !
 * Because the destructor will eventually be called on the objects, the user must initialize them if they are not
 * trivially destructible (eg. containing pointers of any sort).
 */
template<typename Arr>
unique_malloc_ptr<typename concepts::extent_kind<Arr>::unknown_bound_array[]>
make_unique_aligned_uninitialized(size_t n)
{
    using base_type = typename concepts::extent_kind<Arr>::base_type;
    return unique_malloc_ptr<base_type[]>(details::free_delete<base_type[]>::allocate(n));
}

/** A static array (on the heap or stack) of unitialized elements
 *
 * The user is responsable of only acessing elements initialized with emplace_at().
 * Before destruction all initialized elements must be deleted with pop_at()
 */
template<typename T, size_t size, bool heap = true, size_t align = alignof(T)> class uninitialized_array;

namespace details {

template<typename T, size_t size, bool heap, size_t align> class static_array_interface
{
    using Base = uninitialized_array<T, size, heap, align>;

  public:
    template<typename... Args> T& emplace_at(size_t pos, Args&&... args)
    {
        T& ref = this->operator[](pos);
        new (&ref) T(std::forward<Args>(args)...);
        return ref;
    }

    T pop_at(size_t pos)
    {
        T& ref             = this->operator[](pos);
        T              res = std::move(ref);
        ref.~T();
        return res;
    }

    T& operator[](size_t pos)
    {
        check_idx(pos);
        return *reinterpret_cast<T*>(&static_cast<Base*>(this)->_arr[pos]);
    }

    const T& operator[](size_t pos) const
    {
        check_idx(pos);
        return *reinterpret_cast<const T*>(&static_cast<const Base*>(this)->_arr[pos]);
    }

  protected:
    static void check_idx(size_t i) { assume(i < size, "Index out of bound %ull >= %ull", i, size); }
    static_assert(align >= alignof(T), "Alignment requirement must be greater or equal to the contained type");
};

} /* namespace details */

template<typename T, size_t size, size_t align>
class uninitialized_array<T, size, true, align> : public details::static_array_interface<T, size, true, align>
{
    friend class details::static_array_interface<T, size, true, align>;
    using storage_t = typename std::aligned_storage<sizeof(T), align>::type;
    unique_malloc_ptr<storage_t[]> _arr;

  public:
    uninitialized_array()
      : _arr(make_unique_aligned_uninitialized<storage_t[]>(size))
    {}
    uninitialized_array(uninitialized_array&&) = default;
    uninitialized_array& operator=(uninitialized_array&&) = default;
    uninitialized_array(const uninitialized_array&)       = delete;
    uninitialized_array& operator=(const uninitialized_array&) = delete;
};

template<typename T, size_t size, size_t align>
class uninitialized_array<T, size, false, align> : public details::static_array_interface<T, size, false, align>
{
    friend class details::static_array_interface<T, size, false, align>;
    typename std::aligned_storage<sizeof(T), align>::type _arr[size];

  public:
    uninitialized_array(uninitialized_array&&) = default;
    uninitialized_array& operator=(uninitialized_array&&) = default;
    uninitialized_array(const uninitialized_array&)       = delete;
    uninitialized_array& operator=(const uninitialized_array&) = delete;
};

} /* namespace unsafe */

/** Allocates an object of type T, respecting overalignement requirement
 */
template<typename T, typename... Args>
unique_malloc_ptr<typename concepts::extent_kind<T>::single_object>
make_unique_aligned(Args&&... args)
{
    auto ptr = unsafe::make_unique_aligned_uninitialized<T>();
    new (ptr.get()) T(std::forward<Args>(args)...);
    return ptr;
}

/** Allocates n objects of type T, respecting overalignement requirement
 * The object are default initialized.
 */
template<typename Arr>
unique_malloc_ptr<typename concepts::extent_kind<Arr>::unknown_bound_array[]>
make_unique_aligned(size_t n)
{
    using base_type = typename concepts::extent_kind<Arr>::base_type;

    auto uptr = unsafe::make_unique_aligned_uninitialized<Arr>(n);
    new (uptr.get()) base_type[n]{};
    return uptr;
}
}

#endif // MEMORY_HPP
