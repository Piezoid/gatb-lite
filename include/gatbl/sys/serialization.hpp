/// Defines serialization function for contiguous storage types
#ifndef GATBL_RANGE_IO_HPP
#define GATBL_RANGE_IO_HPP

#include <iosfwd>
#include <vector>
#include "gatbl/utils/iterator_pair.hpp"

namespace gatbl {

// Conversions of typed spans to byte spans

using bytes       = span<byte>;
using const_bytes = span<const byte>;

// for pointers

template<typename T>
static inline enable_if_t<std::is_trivially_copyable<T>::value, byte*>
as_bytes(T* p)
{
    return reinterpret_cast<byte*>(p);
}

template<typename T, typename = enable_if_t<std::is_trivially_copyable<T>::value>>
static inline const enable_if_t<std::is_trivially_copyable<T>::value, byte*>
as_bytes(const T* p)
{
    return reinterpret_cast<const byte*>(p);
}

// for refererences

template<typename T>
static inline enable_if_t<std::is_trivially_copyable<T>::value, bytes>
as_bytes(T& ref)
{
    return {reinterpret_cast<byte*>(&ref), sizeof(T)};
}

template<typename T>
static inline enable_if_t<std::is_trivially_copyable<T>::value, const_bytes>
as_bytes(const T& ref)
{
    return {reinterpret_cast<const byte*>(&ref), sizeof(T)};
}

// for spans

template<typename T>
static inline enable_if_t<std::is_trivially_copyable<T>::value, bytes>
as_bytes(span<T> r)
{
    return {reinterpret_cast<byte*>(r.begin()), r.size() * sizeof(T)};
}

template<typename T>
static inline enable_if_t<std::is_trivially_copyable<T>::value, const_bytes>
as_bytes(span<const T> r)
{
    return {reinterpret_cast<const byte*>(r.begin()), r.size() * sizeof(T)};
}

// Read/write for iostream

template<typename T, typename CharT, typename Traits>
static inline auto
write(std::basic_ostream<CharT, Traits>& outs, const T& v)
  -> decltype(concepts::type_require<std::basic_ostream<CharT, Traits>&>(as_bytes(v)))
{
    static_assert(sizeof(T) % sizeof(CharT) == 0, "invalid char type");
    auto span = as_bytes(v);
    outs.write(reinterpret_cast<const CharT*>(span.begin()), std::streamsize(span.size() / sizeof(CharT)));
    return outs;
}

template<typename T, typename CharT, typename Traits>
static inline auto
read(std::basic_istream<CharT, Traits>& ins, T&& v)
  -> decltype(concepts::type_require<std::basic_istream<CharT, Traits>&>(as_bytes(v)))
{
    static_assert(sizeof(T) % sizeof(CharT) == 0, "invalid char type");
    auto span = as_bytes(v);
    ins.read(reinterpret_cast<CharT*>(span.begin()), std::streamsize(span.size() / sizeof(CharT)));
    return ins;
}

// read/write adapter for std::vectors

template<typename O, typename T>
static inline enable_if_t<std::is_trivially_copyable<T>::value, O&>
write(O& out, const std::vector<T>& v)
{
    write(out, v.size());
    write(out, span<const T>(v.data(), v.size()));
    return out;
}

template<typename I, typename T>
static inline enable_if_t<std::is_trivially_copyable<T>::value, I&>
read(I& in, std::vector<T>& v)
{
    size_t size;
    read(in, size);
    v.resize(size);
    read(in, span<T>(v.data(), v.size()));
    return in;
}

} // namespace gatbl

#endif // GATBL_RANGE_IO_HPP
