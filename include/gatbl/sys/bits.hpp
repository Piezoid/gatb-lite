#ifndef GATBL_BITS_HPP
#define GATBL_BITS_HPP

#include "gatbl/common.hpp"
#include <climits>
#include <cstdint>
#include <limits>
#include <utility>
#include <algorithm>

#include "gatbl/utils/ranges.hpp"

namespace gatbl {

namespace concepts {

template<typename T, typename Res_t = T> using if_unsigned_t = enable_if_t<(T(-1) > T(0)), Res_t>;

template<typename T> auto Unsigned(T) -> if_unsigned_t<T>;

/// Checks for a unsigned integer or a finite unsigned random range
template<typename T> auto RawKmer(T) -> if_unsigned_t<T>;
template<typename Range, typename T = Range&&>
auto
RawKmer(Range&& r) -> decltype(type_require<T&&>(Unsigned(at(r, 0u)), size(r)));

/// Checks for a finite unsigned random range
template<typename Range, typename T = Range&&>
auto
UnsignedRange(Range&& r) -> decltype(type_require<T&&>(Unsigned(at(r, 0u)), size(r)));

/// Checks that an unsigned random access range is mutable
template<typename Range, typename T = Range&&>
auto
UnsignedMutableRange(Range&& r) -> decltype(type_require<Range&&>(Unsigned(at(r, 0u) = at(r, 0u)), size(r)));

/// Checks that one unsigned random access range content is assigneable to a second range
template<typename Rin, typename Rout>
auto
CopyableUnsignedRanges(const Rin& in, Rout&& out) -> decltype(value_require(std::forward<Rout>(out),
                                                                            Unsigned(*begin(in)),
                                                                            Unsigned(at(out, 0u) = at(in, 0u)),
                                                                            size(in),
                                                                            size(out)));
} // namespace concepts

namespace bits {

namespace details {
template<size_t bits> struct uint;

template<> struct uint<8>
{
    using type      = uint8_t;
    using fast_type = uint_fast8_t;
};

template<> struct uint<16>
{
    using type      = uint16_t;
    using fast_type = uint_fast16_t;
};

template<> struct uint<32>
{
    using type      = uint32_t;
    using fast_type = uint_fast32_t;
};

template<> struct uint<64>
{
    using type      = uint64_t;
    using fast_type = uint_fast64_t;
};

/// byte order swaping functions
inline constexpr uint16_t
bswap(uint16_t v)
{
    return __builtin_bswap16(v);
}
inline constexpr uint32_t
bswap(uint32_t v)
{
    return __builtin_bswap32(v);
}
inline constexpr uint64_t
bswap(uint64_t v)
{
    return __builtin_bswap64(v);
}

#ifdef __SIZEOF_INT128__
template<> struct uint<128>
{
    using type      = __uint128_t;
    using fast_type = __uint128_t;
};
static constexpr size_t max_uint_bitwidth = 128;

union splited_int128_t
{
    __uint128_t full;
    uint64_t    splited[2];
};

inline CPP14_CONSTEXPR __uint128_t
                       bswap(__uint128_t v)
{ // FIXME:: use suffle
    splited_int128_t s{v};
    uint64_t         tmp = bswap(s.splited[0]);
    s.splited[0]         = bswap(s.splited[1]);
    s.splited[0]         = tmp;
    return s.full;
}

#else
static constexpr size_t max_uint_bitwidth = 64;
#endif

#define BITSET_DEFAULT_128 1

#if defined(__SIZEOF_INT128__) && BITSET_DEFAULT_128
using bitset_default_word_t = __uint128_t;
#else
using bitset_default_word_t               = size_t;
#endif

} // namespace details

template<size_t N> using uint_t      = typename details::uint<N>::type;
template<size_t N> using uint_fast_t = typename details::uint<N>::fast_type;

template<typename T>
inline constexpr concepts::if_unsigned_t<T, std::integral_constant<size_t, sizeof(T) * CHAR_BIT>>
bitwidth()
{
    return {};
}

template<typename T>
inline constexpr concepts::if_unsigned_t<T, std::integral_constant<size_t, sizeof(T) * CHAR_BIT>>
bitwidth(const T&)
{
    return {};
}

template<typename R>
inline constexpr auto
bitwidth(const R& r) -> decltype(bitwidth<decltype(*begin(r))>() * size(r))
{
    return bitwidth<decltype(*begin(r))>() * size(r);
}

template<typename R>
inline constexpr auto
bitwidth() -> decltype(bitwidth<value_t<R>>() * size(R{}))
{
    return bitwidth<value_t<R>>() * size(R{});
}

template<typename T>
constexpr inline concepts::if_unsigned_t<T, std::pair<size_t, size_t>>
split_bitoffset(size_t n)
{
    constexpr size_t word_bits  = bitwidth<T>();
    size_t           word_shift = n / word_bits;
    size_t           bit_shift  = n - (word_shift * word_bits);
    return {word_shift, bit_shift};
}

inline constexpr bool
parity(unsigned int v)
{
    return __builtin_parity(v) == 1;
}

inline constexpr bool
parity(unsigned long v)
{
    return __builtin_parityl(v) == 1;
}

inline constexpr bool
parity(unsigned long long v)
{
    return __builtin_parityll(v) == 1;
}

/// Returns true for odd parity
template<typename R>
inline CPP14_CONSTEXPR auto
parity(const R& r) -> decltype(concepts::type_require<bool>(concepts::UnsignedRange(r)))
{
    bool p = false;
    for (const auto& v : r) {
        p = p != parity(v);
    }
    return p;
}

inline CPP14_CONSTEXPR size_t
                       clz(unsigned int x)
{
    return uint8_t(__builtin_clz(x));
}

inline CPP14_CONSTEXPR size_t
                       clz(unsigned long x)
{
    return uint8_t(__builtin_clzl(x));
}

inline CPP14_CONSTEXPR size_t
                       clz(unsigned long long x)
{
    return uint8_t(__builtin_clzll(x));
}

/// CLZ for uint smaller than 32 bits
template<typename T>
inline CPP14_CONSTEXPR auto
clz(T x) -> enable_if_t<sizeof(T) < sizeof(unsigned int), decltype(clz(static_cast<unsigned int>(x)))>
{

    auto         lz         = clz(static_cast<unsigned int>(x));
    decltype(lz) delta_bits = CHAR_BIT * (sizeof(unsigned int) - sizeof(T));
    assume(lz >= delta_bits, "Less leading zeros than expected");
    return lz - delta_bits;
}

template<typename T>
inline CPP14_CONSTEXPR auto
clz(T x)
  -> enable_if_t<(sizeof(T) > sizeof(unsigned long long)) && sizeof(T) % sizeof(unsigned long long) == 0, unsigned>
{
    constexpr size_t n = sizeof(T) / sizeof(unsigned long long);

    auto   y      = reinterpret_cast<unsigned long long(&)[n]>(x);
    size_t result = 0;
    for (size_t i = n; i-- > 0;) {
        if (y[i] == 0)
            result += sizeof(T) * CHAR_BIT;
        else {
            result += clz(y[i]);
            break;
        }
    }
    return result;
}

/// Returns the log2 rounded up
template<typename T>
inline CPP14_CONSTEXPR auto
ilog2(const T& v) -> decltype(bitwidth<T>() - clz(v - 1))
{
    constexpr T one = 1;
    return v > one ? bitwidth<T>() - clz(v - T(1u)) : 0;
}

/// Returns the log2(x+1) rounded up, this is the number of bits required to store the number
template<typename T>
inline CPP14_CONSTEXPR auto
ilog2p1(const T& v) -> decltype(bitwidth<T>() - clz(v))
{
    constexpr T one = 1;
    return v > one ? bitwidth<T>() - clz(v) : 1;
}

/// Round up to next power of two, greater or equal to the input
template<typename T>
inline CPP14_CONSTEXPR auto
clp2(size_t v) -> decltype(T(ilog2(v)))
{
    return T(1) << ilog2(v);
}

template<typename R>
constexpr auto
rshift(R& r, size_t n) -> decltype(rshift(r, r, n))
{
    return rshift(r, r, n);
}

template<typename Rin, typename Rout>
inline CPP14_CONSTEXPR auto
rshift(const Rin& in, Rout& out, size_t n) -> decltype(concepts::CopyableUnsignedRanges(in, out))
{
    // Describe the range element
    using T                    = decltype(at(in, 0));
    constexpr size_t word_bits = bitwidth<T>();

    const auto bit_offset = split_bitoffset<T>(n);
    const auto word_shift = bit_offset.first;
    const auto bit_shift  = bit_offset.second;
    const auto szi        = size(in);
    const auto szo        = size(out);
    assume(szi <= szo, "rshift: size(in)=%lu > size(out)%lu", szi, szo);
    assume(word_shift < szo, "word_shift=%lu >= size(out)=%lu", word_shift, szi);

    size_t i = 0;
    if (bit_shift != 0) {
        for (; i < szi - word_shift; i++) {
            const T src = at(in, szi - 1 - word_shift - i);
            if (i > 0) { at(out, szi - i) |= src << (word_bits - bit_shift); }
            at(out, szi - 1 - i) = src >> bit_shift;
        }
    } else {
        for (; i < szi - word_shift; i++) {
            at(out, szi - 1 - i) = at(in, szi - 1 - word_shift - i);
        }
    }

    for (; i < szo; i++) {
        at(out, szi - 1 - i) = 0;
    }

    return out;
}

template<typename Rin, typename Rout>
inline CPP14_CONSTEXPR auto
rshift_restricted(const Rin& in, Rout& out, size_t n) -> decltype(concepts::CopyableUnsignedRanges(in, out))
{
    using T                    = decltype(at(in, 0));
    constexpr size_t word_bits = bitwidth<T>();

    const auto bit_offset = split_bitoffset<T>(n);
    const auto word_shift = bit_offset.first;
    const auto bit_shift  = bit_offset.second;
    const auto szi        = size(in);
    const auto szo        = size(out);
    assume(szi <= szo, "size(in)=%lu > size(out)%lu", szi, szo);
    assume(word_shift < szo, "word_shift=%lu >= size(out)=%lu", word_shift, szi);

    size_t i = 0;
    for (; i < word_shift; i++) {
        at(out, i) = T(0u);
    }

    if (bit_shift != 0) {
        T remainder{0u};
        for (; i < szi; i++) {
            const T& src = at(in, i - word_shift);
            at(out, i)   = (src >> bit_shift) | remainder;
            remainder    = src << (word_bits - bit_shift);
        }
    } else {
        for (; i < szi; i++) {
            at(out, i) = at(in, i - word_shift);
        }
    }
    return std::forward<Rout>(out);
}

template<typename R>
constexpr auto
lshift(R&& r, size_t n) -> decltype(lshift(r, std::forward<R>(r), n))
{
    return lshift(r, r, n);
}

template<typename Rin, typename Rout>
CPP14_CONSTEXPR auto
lshift(const Rin& in, Rout& out, size_t n) -> decltype(concepts::CopyableUnsignedRanges(in, out))
{
    using T                    = decltype(at(in, 0));
    constexpr size_t word_bits = bitwidth<T>();

    const auto bit_offset = split_bitoffset<T>(n);
    const auto word_shift = bit_offset.first;
    const auto bit_shift  = bit_offset.second;
    const auto szi        = size(in);
    const auto szo        = size(out);
    assume(szi <= szo, "size(in)=%lu > size(out)%lu", szi, szo);
    assume(word_shift < szo, "word_shift=%lu >= size(out)=%lu", word_shift, szi);

    size_t i = 0;
    if (bit_shift != 0) {
        for (; i < szi - word_shift; i++) {
            const T& src = at(in, i + word_shift);
            if (i > 0) { at(out, i - 1) |= src >> (word_bits - bit_shift); }
            at(out, i) = src << bit_shift;
        }
    } else {
        for (; i < szi - word_shift; i++) {
            at(out, i) = at(in, i + word_shift);
        }
    }

    for (; i < szo; i++) {
        at(out, i) = 0;
    }

    return out;
}

template<typename Rin, typename Rout>
CPP14_CONSTEXPR auto
lshift_restricted(const Rin& in, Rout& out, size_t n) -> decltype(concepts::CopyableUnsignedRanges(in, out))
{
    using T                    = decltype(at(in, 0));
    constexpr size_t word_bits = bitwidth<T>();

    const auto bit_offset = split_bitoffset<T>(n);
    const auto word_shift = bit_offset.first;
    const auto bit_shift  = bit_offset.second;
    const auto szi        = size(in);
    const auto szo        = size(out);
    assume(szi <= szo, "size(in)=%lu > size(out)%lu", szi, szo);
    assume(word_shift < szo, "word_shift=%lu >= size(out)=%lu", word_shift, szi);

    size_t i = 0;
    for (; i < word_shift; i++) {
        at(out, szi - 1 - i) = T(0u);
    }

    if (bit_shift != 0) {
        T remainder{0u};
        for (; i < szi; i++) {
            const T& src         = at(in, szi - 1 + word_shift - i);
            at(out, szi - 1 - i) = (src << bit_shift) | remainder;
            remainder            = src >> (word_bits - bit_shift);
        }
    } else {
        for (; i < szi; i++) {
            at(out, szi - 1 - i) = at(in, szi - 1 + word_shift - i);
        }
    }

    return out;
}

template<typename T>
inline constexpr concepts::if_unsigned_t<T>
rotl(T x, size_t n)
{
    constexpr size_t bits = bitwidth<T>();
    assume(n < bits, "shift=%lu >= bitwidth=%lu", n, bits);
    return (x << n) | (x >> (bits - n));
}

template<typename T>
inline constexpr concepts::if_unsigned_t<T>
rotr(T x, size_t n)
{
    constexpr size_t bits = bitwidth<T>();
    assume(n < bits, "shift=%lu >= bitwidth=%lu", n, bits);
    return (x >> n) | (x << (bits - n));
}

/// Resilient bitmask that accept parameters producing clipped bitmask
template<typename T>
inline CPP14_CONSTEXPR concepts::if_unsigned_t<T>
                       bitmask(size_t n, size_t pos = 0)
{
    constexpr size_t bits = bitwidth<T>();
    assume(n <= bits, "n=%lu > bitwidth=%lu", n, bitwidth<T>());
    assume(pos + n <= bits, "pos=%lu + n=%lu > bitwidth=%lu", pos, n, bitwidth<T>());
    return (T(-1) >> (bits - n)) << pos;
}

/// A bitmask assuming parameters for unclipped output
template<typename T>
inline CPP14_CONSTEXPR concepts::if_unsigned_t<T>
                       bitmask_unclipped(size_t n, size_t pos = 0)
{
    constexpr size_t bits = bitwidth<T>();
    assume(n > 0, "zero weight mask");
    assume(n <= bits, "mask weight's %lu > capacity=%lu", n, bits);
    assume(pos < bits, "pos=%lu >= capacity%lu", pos, bits);
    assume(n + pos <= bits, "pos+weight=%lu > capacity=%lu", n + pos, bits);
    return (T(-1) >> (bits - n)) << pos;
}

template<typename T>
inline constexpr concepts::if_unsigned_t<T>
peek_bitgroup(T on, size_t n = 1, size_t pos = 0)
{
    return (on >> pos) & bitmask<T>(n);
}

template<typename Range, typename Result = value_t<Range>>
inline auto
peek_bitgroup(const Range& r, size_t n = 1, size_t pos = 0) noexcept -> decltype(concepts::value_require<Result>(
  value_t<Range>{at(r, size(r) - 1)},
  enable_if_t<(sizeof(Result) <= sizeof(at(r, 0))) && (Result(-1) > Result(0)), Result>{}))
{
    using range_element_t             = value_t<Range>;
    static constexpr size_t word_bits = bits::bitwidth<range_element_t>();
    static constexpr size_t T_nbits   = bits::bitwidth<Result>();
    assume(n <= T_nbits, "mask weight=%lu > word capacity=%lu");

    const auto   bit_offset      = split_bitoffset<range_element_t>(pos);
    const auto   word_shift      = bit_offset.first;
    const auto   bit_shift       = bit_offset.second;
    const size_t shifted_end_bit = bit_shift + n;
    assume(word_shift + (shifted_end_bit > word_bits) < size(r),
           "peek bits past the end at word %lu >= size(r)=%lu",
           word_shift + (shifted_end_bit > word_bits),
           size(r));

    if (shifted_end_bit <= word_bits) {
        return peek_bitgroup(at(r, word_shift), n, bit_shift);
    } else {
        const size_t lowbits  = word_bits - bit_shift;
        const size_t highbits = shifted_end_bit - word_bits;
        Result       res      = peek_bitgroup(at(r, word_shift), lowbits, bit_shift);
        res |= Result{peek_bitgroup(at(r, word_shift + 1), highbits, 0)} << lowbits;
        return res;
    }
}

/// Set bits pos...pos+n-1 in `target` with bits 0..n-1 from `from`
/// If all bits are 0 in the target region, a bitwise OR would be faster, even more so if `from` is already masked.
template<typename T>
inline constexpr concepts::if_unsigned_t<T>&
poke_bitgroup(T& target, T from, size_t n = 1, size_t pos = 0)
{
    T mask = bitmask<T>(n, pos);
    target ^= (target ^ (from << pos)) & mask;
    return target;
}

template<typename Range>
inline constexpr auto
poke_bitgroup(Range& r, value_t<Range> v, size_t n = 1, size_t pos = 0) -> decltype(concepts::UnsignedRange(r))
{
    using range_element_t      = value_t<Range>;
    constexpr size_t word_bits = bits::bitwidth<range_element_t>();
    assume(n <= word_bits, "mask weight=%lu > word capacity=%lu");

    const auto   bit_offset      = split_bitoffset<range_element_t>(pos);
    const auto   word_shift      = bit_offset.first;
    const auto   bit_shift       = bit_offset.second;
    const size_t shifted_end_bit = bit_shift + n;
    assume(word_shift + (shifted_end_bit > word_bits) < size(r),
           "peek bits past the end at word %lu >= size(r)=%lu",
           word_shift + (shifted_end_bit > word_bits),
           size(r));

    if (shifted_end_bit <= word_bits) {
        at(r, word_shift) = poke_bitgroup(at(r, word_shift), range_element_t(v), n, bit_shift);
    } else {
        const size_t lowbits  = word_bits - bit_shift;
        const size_t highbits = shifted_end_bit - word_bits;
        at(r, word_shift)     = poke_bitgroup(at(r, word_shift), range_element_t(v), lowbits, bit_shift);
        at(r, word_shift + 1) = bits::poke_bitgroup(at(r, word_shift + 1), range_element_t(v >> lowbits), highbits, 0);
    }

    return r;
}

/// Reverse order of bit groups in log_2(N) where N is the number of group per word (must be a power of two)
template<size_t group_bits = 1, typename T>
inline CPP14_CONSTEXPR concepts::if_unsigned_t<T>
                       reverse_bitgroups(T v)
{
    static_assert((group_bits & (group_bits - 1)) == 0, "The group size must be a power of two");
    size_t level = bitwidth<T>();
    T      mask  = ~T(0);

    // First, try to use specialized instr for byte order reversal
    CPP17_IF_CONSTEXPR(group_bits <= 8 && bitwidth<T>() > 1)
    {
        v = details::bswap(v);
        // Update the mask for byte width: 0xff00ff00...
        while (level > CHAR_BIT) {
            level >>= 1;
            mask ^= (mask << level);
        }
    }

    while (level > group_bits) {
        level >>= 1;
        mask ^= (mask << level);
        v = ((v >> level) & mask) | ((v << level) & ~mask);
    }

    return v;
}

template<size_t group_bits = 1, typename R>
inline CPP14_CONSTEXPR auto
reverse_bitgroups(R&& r) -> decltype(concepts::UnsignedMutableRange(r))
{
    std::reverse(begin(r), end(r));
    for (auto& w : r) {
        w = reverse_bitgroups<group_bits>(w);
    }
    return std::forward<R>(r);
}

template<size_t group_bits = 1, typename Rin, typename Rout>
inline CPP14_CONSTEXPR auto
reverse_bitgroups(const Rin& in, Rout&& out) -> decltype(concepts::CopyableUnsignedRanges(in, std::forward<Rout>(out)))
{
    assume(size(in) <= size(out), "size(in)=%lu > size(out)=%lu", size(in), size(out));
    std::copy(rbegin(in), rend(in), begin(out));
    for (auto& w : out) {
        w = reverse_bitgroups<group_bits>(w);
    }
    return std::forward<Rout>(out);
}

template<size_t group_bits = 1, typename T>
inline CPP14_CONSTEXPR auto
tile_bitgroup(T v) -> concepts::if_unsigned_t<T>
{
    constexpr size_t bits  = bitwidth<T>();
    size_t           width = group_bits;
    while (width < bits) {
        v |= v << width;
        width <<= 1;
    }
    return v;
}

} // namespace bits

} // namespace gatbl

#endif // GATBL_BITS_HPP
