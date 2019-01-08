#ifndef MMAP_HPP
#define MMAP_HPP

#include <unistd.h>
#include <sys/mman.h>
#include "gatbl/common.hpp"
#include "gatbl/utils/unique_range.hpp"
#include "gatbl/sys/exceptions.hpp"

namespace gatbl {

namespace details {
// FIXME: < c++17
static inline const size_t page_size = sys::check_ret(::sysconf(_SC_PAGESIZE), "page size");

template<typename T> struct munmapper
{
    void operator()(T* addr, size_t n) const
    {
        const size_t bytes = size_to_bytes(n);
        if (unlikely(bytes == 0)) return;
        sys::check_ret(::munmap(unconst_void(addr), bytes), "munmap");
    }

    static constexpr inline void* unconst_void(const void* p) noexcept { return const_cast<void*>(p); }

    static inline size_t size_to_bytes(size_t n)
    {
        const size_t bytes = sizeof(T) * n;
        // Round up to page size
        return ((bytes + page_size - 1) / page_size) * page_size;
    }
};

} // namespace detail

template<typename T> class mmap_range : public unique_range<T, details::munmapper<T>>
{
    using unmapper = details::munmapper<T>;
    using base     = unique_range<T, unmapper>;
    using base::base;

  public:
    using base::operator=;
    mmap_range(int                     fd,
               size_t                  len,
               int                     prot   = PROT_READ,
               int                     flags  = MAP_PRIVATE,
               std::remove_const_t<T>* addr   = nullptr,
               off_t                   offset = 0)
      : base()
    {
        const size_t bytes = details::munmapper<T>::size_to_bytes(len);
        T*           ptr   = reinterpret_cast<T*>(::mmap(addr, bytes, prot, flags, fd, offset));
        if (likely(ptr != MAP_FAILED)) {
            *this = base(ptr, len, {});
        } else {
            sys::throw_syserr("mmap(%p, %lu, %d, %d, %d, %ld)", addr, bytes, prot, flags, fd, offset);
        }
    }

    void advise_sequential(size_t len = 0, size_t from = 0)
    {
        this->advise(MADV_SEQUENTIAL, len, from, "madvise(MADV_SEQUENTIAL, %2$lu, %3$d)");
    }

    void advise_hugepage(size_t len = 0, size_t from = 0)
    {
        this->advise(MADV_HUGEPAGE, len, from, "madvise(MADV_HUGE_PAGE, %2$lu, %3$d)");
    }

  protected:
    void advise(int advice, size_t len = 0, size_t from = 0, const char error_msg[] = "madvise(%p, %lu, %d)")
    {
        len = len > 0 ? len : this->size();
        assume(
          from + len <= this->size(), "advise past the end (offset+len=%lu, size()=%lu)", from + len, this->size());

        void*  addr  = unmapper::unconst_void(this->begin() + sizeof(T) * from);
        size_t bytes = sizeof(T) * len;

        sys::check_ret(::madvise(addr, bytes, advice), error_msg, addr, bytes, advice);
    }
};

} // namespace gatbl

#endif // MMAP_HPP
