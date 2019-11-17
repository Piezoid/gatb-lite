#ifndef FILE_HPP
#define FILE_HPP

#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "gatbl/utils/compatibility.hpp"
#include "gatbl/sys/serialization.hpp"
#include "gatbl/sys/exceptions.hpp"
#include "gatbl/sys/mmap.hpp"

namespace gatbl {

/// Wrap a POSIX file descriptor into a C++ interface
/// The offset state is managed in userspace
struct file_descriptor
{
    file_descriptor(const std::string& path, int flags = O_RDONLY, mode_t mode = S_IRUSR | S_IWUSR)
      : file_descriptor(check_ret(::open(path.c_str(), flags, mode), "open(\"%s\", %o)", path.c_str(), mode))
    {}

    file_descriptor(file_descriptor&& from) noexcept
      : _fd(from._fd)
      , _size(from._size)
      , _pos(from._pos)
    {
        from._fd = -1;
    }

    /// Move constructor
    file_descriptor& operator=(file_descriptor&& from) noexcept
    {
        close();

        _fd   = from._fd;
        _size = from._size;
        _pos  = from._pos;

        from._fd = -1;

        return *this;
    }

    /// Copy constructor
    file_descriptor(const file_descriptor& from)
      : _fd(check_ret(::dup(from._fd), "dup(%d)", _fd))
      , _size(from._size)
      , _pos(from._pos)
    {}

    /// Assignement constructor
    file_descriptor& operator=(const file_descriptor& from)
    {
        close();

        _fd   = check_ret(dup(from._fd), "dup(%d)", from._fd);
        _size = from._size;
        _pos  = from._pos;

        return *this;
    }

    static file_descriptor make_shm(const char* name   = "tmp",
                                    int         oflag  = O_RDWR | O_CREAT,
                                    mode_t      mode   = 0700,
                                    bool        unlink = true)
    {
        const int fd = check_ret(::shm_open(name, oflag, mode), "shm_open(\"%s\", %o, %o)", name, oflag, mode);
        if (unlink) { check_ret(::shm_unlink(name), "shm_unlink(\"%s\")", name); }
        return file_descriptor(fd);
    }

    static file_descriptor make_memfd(const char* name = "", unsigned int flags = 0)
    {
        const int fd = check_ret(::memfd_create(name, flags), "memfd_create");
        return file_descriptor(fd);
    }

    size_t truncate(size_t length)
    {
        check_ret(::ftruncate(_fd, off_t(length)), "ftruncate(%d, %zd)", _fd, length);
        _size = length;
        return length;
    }

    void advise(int advise, off_t offset = 0, off_t len = 0)
    {
        if (offset == 0 && len == 0) { len = off_t(this->size()); }

#ifdef POSIX_FADV_SEQUENTIAL
        check_ret(
          ::posix_fadvise(_fd, offset, len, advise), "posix_fadvise(%d, %zd, %zd, %d)", _fd, offset, len, advise);
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

    /// Read file content from offset into the buffer
    /// The buffer begin() is advanced by the number of bytes read
    ssize_t pread(span<byte>& buf_bytes, size_t offset)
    {
        ssize_t sz = check_ret(::pread(_fd, buf_bytes.begin(), buf_bytes.size(), off_t(offset)),
                               "pread(%d, buf, %zd, %zd)",
                               _fd,
                               buf_bytes.size(),
                               offset);
        buf_bytes  = {buf_bytes.begin() + sz, buf_bytes.end()};
        return sz;
    }

    /// Write buffer at offset
    /// The buffer begin() is advanced by the number of bytes written
    ssize_t pwrite(span<const byte>& buf_bytes, size_t offset)
    {
        ssize_t sz    = check_ret(::pwrite(_fd, buf_bytes.begin(), buf_bytes.size(), off_t(offset)),
                               "pwrite(%d, buf, %zd, %zd)",
                               _fd,
                               buf_bytes.size(),
                               offset);
        buf_bytes     = {buf_bytes.begin() + sz, buf_bytes.end()};
        size_t endpos = offset + size_t(sz);
        if (endpos > this->_size) _size = endpos;
        return sz;
    }

    /// Read file content from offset into the buffer
    /// The buffer begin() is advanced by the number of bytes read
    ssize_t read(span<byte>& buf)
    {
        ssize_t sz = file_descriptor::pread(buf, this->tell());
        _pos += size_t(sz);
        return sz;
    }

    /// Write buffer at current position
    /// The buffer begin() is advanced by the number of bytes written
    ssize_t write(span<const byte>& buf)
    {
        ssize_t sz = file_descriptor::pwrite(buf, this->tell());
        _pos += size_t(sz);
        return sz;
    }

#ifdef __linux__
    ssize_t readahead(off64_t offset, size_t count) const
    {
        return check_ret(::readahead(_fd, offset, count), "readahead(%d, %zd, %zu)", _fd, offset, count);
    }
#else
    size_t readahead(off_t offset, size_t count) const { return 0; }
#endif

    void close()
    {
        if (_fd >= 0) { check_ret(::close(_fd), "close(%d)", _fd); }
        _fd = -1;
    }

    ~file_descriptor() { close(); }

    operator bool() const noexcept { return _fd > 0; }

    size_t size() const noexcept { return _size; }
    size_t tell() const noexcept { return _pos; }
    bool   eof() const noexcept { return tell() >= size(); }

    /// We use a user-space position since async read don't use the native file descriptor position
    size_t seek(off_t offset, int whence = SEEK_SET)
    {
        switch (whence) {
            case SEEK_SET: break;
            case SEEK_CUR: offset += tell(); break;
            case SEEK_END: offset += size(); break;
        }

        // clamp
        size_t new_pos = likely(offset > 0) ? size_t(offset) : 0;
        new_pos        = likely(new_pos < _size) ? new_pos : _size;
        return _pos    = new_pos;
    }

  protected:
    explicit file_descriptor(int fd)
      : _fd(fd)
    {
        struct stat _stat = {};
        check_ret(::fstat(fd, &_stat), "stat(%d)", fd);
        _size = size_t(_stat.st_size);
    }
    int getfd() const noexcept { return _fd; }

  protected:
    int    _fd   = -1;
    size_t _size = 0;
    size_t _pos  = 0;
};

template<typename T>
static inline auto
write(file_descriptor& fd, const T& v) -> decltype(concepts::type_require<file_descriptor&>(as_bytes(v)))
{
    const_bytes span = as_bytes(v);
    do
        fd.write(span);
    while (unlikely(!span.empty()));

    return fd;
}

template<typename T>
static inline auto
read(file_descriptor& fd, T&& v) -> decltype(concepts::type_require<file_descriptor&>(as_bytes(v)))
{
    bytes span = as_bytes(v);
    do {
        ssize_t sz = fd.read(span);
        if (sz == 0) throw_syserr("read(): unexpected eof");
    } while (unlikely(!span.empty()));

    return fd;
}

} // namespace gatbl

#endif // FILE_HPP
