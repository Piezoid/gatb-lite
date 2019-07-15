#ifndef FILE_HPP
#define FILE_HPP

#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "gatbl/utils/compatibility.hpp"
#include "gatbl/utils/iterator_pair.hpp"
#include "gatbl/sys/exceptions.hpp"
#include "gatbl/sys/mmap.hpp"

// #include <iostream>
// #include <string>
// #include <vector>
// #include <deque>

// #include "gatbl/ext/filesystem.hpp"
// #include "gatbl/common.hpp"
// #include "gatbl/memory.hpp"
// #include "gatbl/threads.hpp"
//#include<list>
// std::list<int>::iterator

namespace gatbl { namespace sys {

struct sized
{
    size_t size() const noexcept { return _size; }

  protected:
    void setSize(size_t size) noexcept { _size = size; }

  private:
    size_t _size = 0;
};

struct cursor
{
    size_t tell() const noexcept { return _pos; }

  protected:
    size_t setPosition(size_t pos) noexcept { return _pos = pos; }
    size_t incPosition(ssize_t delta) noexcept { return setPosition(tell() + delta); }
    void   operator++() noexcept { incPosition(1); }

  private:
    size_t _pos = 0;
};

struct bound_cursor
  : public sized
  , public cursor
{
    bool eof() const noexcept { return tell() >= size(); }

    /// We use a user-space position since async read don't use the native file descriptor position
    size_t seek(off_t offset, int whence = SEEK_SET)
    {
        switch (whence) {
            case SEEK_SET: break;
            case SEEK_CUR: offset += tell(); break;
            case SEEK_END: offset += size(); break;
        }

        return this->setPosition(offset);
    }

  protected:
    size_t setPosition(ssize_t pos) noexcept
    {
        // clamp
        size_t new_pos = likely(pos > 0) ? size_t(pos) : 0;
        new_pos        = likely(new_pos < size()) ? new_pos : size();
        return cursor::setPosition(new_pos);
    }
    size_t incPosition(ssize_t delta) noexcept { return setPosition(ssize_t(tell()) + delta); }
    void   operator++() noexcept { incPosition(1); }
};

struct file_descriptor : public bound_cursor
{
    file_descriptor(const std::string& path, int flags = O_RDONLY, mode_t mode = S_IRUSR | S_IWUSR)
      : file_descriptor(sys::check_ret(::open(path.c_str(), flags, mode), "open"))
    {}

    file_descriptor(file_descriptor&& from) noexcept
      : bound_cursor(std::move(from))
      , _stat(from._stat)
      , _fd(from._fd)
    {
        from._fd = -1;
    }

    /// Move constructor
    file_descriptor& operator=(file_descriptor&& from) noexcept
    {
        close();

        _fd      = from._fd;
        _stat    = from._stat;
        from._fd = -1;

        bound_cursor::operator=(std::move(from));
        return *this;
    }

    /// Copy constructor
    file_descriptor(const file_descriptor& from)
      : bound_cursor(from)

      , _stat(from._stat)
      , _fd(int(sys::check_ret(::dup(from._fd), "dup")))
    {}

    /// Assignement constructor
    file_descriptor& operator=(const file_descriptor& from)
    {
        close();

        _fd   = int(sys::check_ret(dup(from._fd), "dup"));
        _stat = from._stat;

        bound_cursor::operator=(from);
        return *this;
    }

    static file_descriptor make_shm(const char* name   = "tmp",
                                    int         oflag  = O_RDWR | O_CREAT,
                                    mode_t      mode   = 0700,
                                    bool        unlink = true)
    {
        const unsigned fd = sys::check_ret(::shm_open(name, oflag, mode), "shm_open");
        if (unlink) { sys::check_ret(::shm_unlink(name), "shm_unlink"); }
        return file_descriptor(fd);
    }

    size_t blksize() const noexcept { return static_cast<size_t>(_stat.st_blksize); }

    size_t truncate(size_t length)
    {
        sys::check_ret(::ftruncate(_fd, off_t(length)), "ftruncate");
        this->setSize(length);
        return length;
    }

    void advise(int advise, off_t offset = 0, off_t len = 0)
    {
        if (offset == 0 && len == 0) { len = off_t(this->size()); }

#ifdef POSIX_FADV_SEQUENTIAL
        sys::check_ret(::posix_fadvise(_fd, offset, len, advise), "posix_fadvise");
#endif
    }

    template<typename T = byte>
    mmap_range<T> mmap(off_t              offset = 0,
                       size_t             len    = 0,
                       int                prot   = PROT_READ,
                       int                flags  = MAP_PRIVATE,
                       remove_const_t<T>* addr   = nullptr)
    {
        if (offset == 0 && len == 0) { len = this->size() / sizeof(T); }
        return mmap_range<T>(_fd, len, prot, flags, addr, offset);
    }

    template<typename Buffer>
    auto pread(Buffer& buf, off_t offset) -> decltype(concepts::type_require<size_t>(as_bytes(buf)))
    {
        using gatbl::as_bytes;
        using gatbl::size;
        const auto buf_bytes = as_bytes(buf);
        return sys::check_ret(::pread(_fd, begin(buf_bytes), size(buf_bytes), offset), "pread");
    }

    template<typename Buffer>
    auto pwrite(const Buffer& buf, off_t offset) -> decltype(concepts::type_require<size_t>(as_bytes(buf)))
    {
        using gatbl::as_bytes;
        using gatbl::size;
        const auto buf_bytes = as_bytes(buf);
        size_t     sz        = sys::check_ret(::pwrite(_fd, begin(buf_bytes), size(buf_bytes), offset), "pwrite");
        if (offset + sz > this->size()) setSize(offset + sz);
        return sz;
    }

    template<typename Buffer> size_t read(Buffer& buf)
    {
        size_t sz = file_descriptor::pread(buf, this->tell());
        this->incPosition(ssize_t(sz));
        return sz;
    }

    template<typename Buffer> size_t write(const Buffer& buf)
    {
        size_t sz = file_descriptor::pwrite(buf, this->tell());
        this->incPosition(ssize_t(sz));
        return sz;
    }

#ifdef __linux__
    size_t readahead(off64_t offset, size_t count) const
    {
        return sys::check_ret(::readahead(_fd, offset, count), "readahead");
    }
#else
     size_t readahead(off_t offset, size_t count) const
     {
        return 0;
     }
#endif
    void close()
    {
        if (_fd >= 0) { sys::check_ret(::close(_fd), "close"); }
        _fd = -1;
    }

    ~file_descriptor() { close(); }

    operator bool() const noexcept { return _fd > 0; }

  protected:
    explicit file_descriptor(unsigned fd)
      : _fd(fd)
    {
        sys::check_ret(::fstat(int(fd), &_stat), "stat");
        this->setSize(size_t(_stat.st_size));
    }

    int getfd() const noexcept { return _fd; }

  private:
    struct stat _stat = {};
    int         _fd   = -1;
};

} // namespace sys
} // namespace gatbl

#endif // FILE_HPP
