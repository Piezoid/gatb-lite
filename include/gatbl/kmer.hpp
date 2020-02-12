#ifndef KMER_HPP
#define KMER_HPP

#include <cstdint>
#include <climits>
#include <cmath>
#include <tmmintrin.h>

#include <limits>
#include <string>
#include <stdexcept>

#include "gatbl/common.hpp"
#include "gatbl/utils/nucleotide.hpp"
#include "gatbl/utils/ranges.hpp"
#include "gatbl/utils/no_conversion.hpp"
#include "gatbl/sys/bits.hpp"
#include "gatbl/utils/range_view.hpp"
#include "gatbl/utils/reverse_range.hpp"

namespace gatbl {

using ksize_t     = bits::bitsize_t; // Nucleotide incides inside kmers
using kmer_t      = uint64_t;        // k-mers, uint64_t for k<32, uint64_t for k<64
using minimizer_t = uint32_t;        // Minimizers

struct complement_functor
{
    template<typename T> auto operator()(T n) const -> decltype(complement(n)) { return complement(n); }
};

template<typename R> using reverse_complement_view = range_view<reverse_range<R>, complement_functor>;

template<typename R>
auto
reverse_complement(R&& r)
  -> decltype(concepts::type_require<reverse_complement_view<R>>(concepts::Range(reverse_complement_view<R>(r))))
{

    return reverse_complement_view<R>(r);
}

template<typename R>
auto
reversed(R&& r) -> decltype(concepts::type_require<reverse_range<R>>(concepts::Range(reverse_range<R>(r))))
{
    return reverse_range<R>(r);
}

// Represents the cardinality of a pow2 sized set. Allows div/mod arithmetic operations on indexes.
template<typename T> struct Pow2
{
    using bitsize_t = bits::bitsize_t;
    Pow2(bitsize_t bits)
      : _bits(bits)
    {
        assume(bits < CHAR_BIT * sizeof(T), "Pow2(%u > %u)", unsigned(bits), unsigned(CHAR_BIT * sizeof(T)));
    }

    bitsize_t bits() const { return _bits; }
    T         value() const { return T(1) << _bits; }
    explicit  operator T() const { return value(); }
    T         max() const { return value() - T(1); }

    friend T  operator*(const T& x, const Pow2& y) { return x << y._bits; }
    friend T& operator*=(T& x, const Pow2& y) { return x <<= y._bits; }
    friend T  operator/(const T& x, const Pow2& y) { return x >> y._bits; }
    friend T& operator/=(T& x, const Pow2& y) { return x >>= y._bits; }
    friend T  operator%(const T& x, const Pow2& y) { return x & y.max(); }
    friend T& operator%=(T& x, const Pow2& y) { return x &= y.max(); }
    Pow2&     operator>>=(bitsize_t d)
    {
        _bits -= d;
        return *this;
    }
    Pow2& operator<<=(bitsize_t d)
    {
        _bits += d;
        return *this;
    }
    friend bool operator<(const T& x, const Pow2& y) { return x < y.value(); }
    friend bool operator<=(const T& x, const Pow2& y) { return x < y.value(); }
    friend T    operator+(const T& x, const Pow2& y) { return x + y.value(); }
    friend T&   operator+=(T& x, const Pow2& y) { return x += y.value(); }
    friend T    operator-(const T& x, const Pow2& y) { return x - y.value(); }
    friend T&   operator-=(T& x, const Pow2& y) { return x -= y.value(); }

  private:
    bitsize_t _bits;
};

struct invalid_dna_exception : public std::domain_error
{
    invalid_dna_exception()
      : std::domain_error("Invalid char in DNA")
    {}
};

inline nuc_t hot_fun
nuc2int(char c) // FIXME: rename
{
    //    if((c == 'a' || c == 'c' || c == 't' || c == 'g'
    //           || c == 'A' || c == 'C' || c == 'T' || c == 'G')) {
    //        return (c >> 1) & 3;
    //    }
    nucint_t d = c - 65;
    if (d <= 51) {
        if (0x0008004500080045ull & (1ull << d)) { return as_nuc(nucint_t(c >> 1) & 3u); }
    }
    return nuc_t::N;
}

inline std::string pure_fun hot_fun
kmer2str(kmer_t num, ksize_t k)
{
    std::string  res(k, 0);
    Pow2<kmer_t> anc(2 * (k - 1));
    for (uint i = 0; i < k; i++) {
        uint nuc = num / anc;
        num %= anc;
        anc >>= 2;

        assume(nuc <= 3, "nuc=%u > 3", nuc);
        res[i] = "ACTG"[nuc];
    }
    return res;
}

inline kmer_t pure_fun hot_fun
str2num(const std::string& str)
{
    kmer_t res(0);
    for (uint i = 0; i < str.size(); i++) {
        res <<= 2;
        res |= as_integer(nuc2int(str[i]));
    }
    return res;
}

struct ReversibleHash
{
    inline int32_t operator()(uint32_t x) const pure_fun hot_fun
    {
        x = ((x >> 16) ^ x) * uint32_t(0x2c1b3c6d);
        x = ((x >> 16) ^ x) * uint32_t(0x297a2d39);
        x = ((x >> 16) ^ x);
        return x;
    }

    inline int64_t operator()(uint64_t x) const pure_fun hot_fun
    {
        x = ((x >> 32) ^ x) * uint64_t(0xD6E8FEB86659FD93);
        x = ((x >> 32) ^ x) * uint64_t(0xD6E8FEB86659FD93);
        x = ((x >> 32) ^ x);
        return x;
    }
};

struct ReversibleHashOp
{
    inline int32_t operator()(uint32_t x) const pure_fun hot_fun
    {
        x = ((x >> 16) ^ x) * uint32_t(0x0cf0b109); // PowerMod[0x297a2d39, -1, 2^32]
        x = ((x >> 16) ^ x) * uint32_t(0x64ea2d65);
        x = ((x >> 16) ^ x);
        return x;
    }

    inline int64_t operator()(uint64_t x) pure_fun hot_fun
    {
        x = ((x >> 32) ^ x) * uint64_t(0xCFEE444D8B59A89B);
        x = ((x >> 32) ^ x) * uint64_t(0xCFEE444D8B59A89B);
        x = ((x >> 32) ^ x);
        return x;
    }
};

// It's quite complex to bitshift mmx register without an immediate (constant) count
// See: https://stackoverflow.com/questions/34478328/the-best-way-to-shift-a-m128i
inline __m128i pure_fun hot_fun
mm_bitshift_left(__m128i x, unsigned count)
{
    assume(count < 128, "count=%u >= 128", count);
    __m128i carry = _mm_slli_si128(x, 8);
    if (count >= 64)                              // TODO: bench: Might be faster to skip this fast-path branch
        return _mm_slli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
    // else
    carry = _mm_srli_epi64(carry, 64 - count);

    x = _mm_slli_epi64(x, count);
    return _mm_or_si128(x, carry);
}

inline __m128i pure_fun hot_fun
mm_bitshift_right(__m128i x, unsigned count)
{
    assume(count < 128, "count=%u >= 128", count);
    __m128i carry = _mm_srli_si128(x, 8);
    if (count >= 64) return _mm_srli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
    // else
    carry = _mm_slli_epi64(carry, 64 - count);

    x = _mm_srli_epi64(x, count);
    return _mm_or_si128(x, carry);
}

inline __uint128_t pure_fun hot_fun
rcb(const __uint128_t& in, uint n)
{
    assume(n <= 64, "n=%u > 64", n);
    union kmer_u
    {
        __uint128_t k;
        __m128i     m128i;
        uint64_t    u64[2];
        uint8_t     u8[16];
    };
    kmer_u res = {.k = in};
    static_assert(sizeof(res) == sizeof(__uint128_t), "kmer sizeof mismatch");

    // Swap byte order
    kmer_u shuffidxs = {.u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};
    res.m128i        = _mm_shuffle_epi8(res.m128i, shuffidxs.m128i);

    // Swap nuc order in bytes
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
    const uint64_t c2 = 0x3333333333333333;
    for (uint64_t& x : res.u64) {
        x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
        x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc
        x ^= 0xaaaaaaaaaaaaaaaa;                      // Complement;
    }

    // Realign to the right
    res.m128i = mm_bitshift_right(res.m128i, 128 - 2 * n);
    return res.k;
}

inline uint64_t pure_fun hot_fun
rcb(uint64_t in, uint n)
{
    assume(n <= 32, "n=%u > 32", n);
    // Complement, swap byte order
    uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
    // Swap nuc order in bytes
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
    const uint64_t c2 = 0x3333333333333333;
    res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
    res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

    // Realign to the right
    res >>= 64 - 2 * n;
    return res;
}

inline uint32_t pure_fun hot_fun
rcb(uint32_t in, uint n)
{
    assume(n <= 16, "n=%u > 16", n);
    // Complement, swap byte order
    uint32_t res = __builtin_bswap32(in ^ 0xaaaaaaaa);

    // Swap nuc order in bytes
    const uint32_t c1 = 0x0f0f0f0f;
    const uint32_t c2 = 0x33333333;
    res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
    res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

    // Realign to the right
    res >>= 32 - 2 * n;
    return res;
}

template<typename T>
constexpr inline auto
get_kmer(T kmer) -> decltype(concepts::RawKmer(kmer))
{
    return kmer;
}

template<typename T> using get_kmer_t = decltype(get_kmer(std::declval<T>()));

/// A wrapper for arguments that disable implicit conversions but allowing to extract kmer from decorated kmers (eg.
/// sized_kmer)
template<typename T> struct no_conversion_kmer : public no_conversion<T>
{
    using no_conversion<T>::no_conversion;

    template<typename U, typename = if_is_same_cvref<T, get_kmer_t<U>>>
    no_conversion_kmer(U&& v) // Implicit constructor with wildcard type argument for catching all conversion attempts
      : no_conversion<T>(get_kmer(std::forward<U>(v)))
    {}
};

template<typename kmer_t = kmer_t, typename ksize_t = ksize_t> struct packed_layout sized_kmer;

template<typename KmerT = size_t, typename KSizeT = ksize_t> struct packed_layout kmer_model
{
    using kmer_t  = KmerT;
    using ksize_t = KSizeT;

    kmer_model() = default;
    kmer_model(int k)
      : _k(ksize_t(k))
    {
        assume(k >= 0 && k <= std::numeric_limits<ksize_t>::max(), "k out of range");
    }
                           operator ksize_t() const { return _k; }
    constexpr ksize_t      size() const { return _k; }
    constexpr ksize_t      bitwidth() const { return 2 * _k; }
    CPP14_CONSTEXPR kmer_t cardinality() const
    {
        assume(bitwidth() < bits::bitwidth<kmer_t>(), "cadinality is not representable in kmer type");
        return kmer_t(1) << bitwidth();
    }
    CPP14_CONSTEXPR kmer_t mask() const { return bits::bitmask<kmer_t>(bitwidth()); }
    CPP14_CONSTEXPR kmer_t mask(no_conversion_kmer<kmer_t> kmer) const { return kmer & mask(); }
    // void assume_masked(no_conversion_kmer<kmer_t> kmer) const { assume(kmer <= mask(), "kmer is not masked"); }
    void assert_masked(no_conversion_kmer<kmer_t> kmer) const { assert(kmer <= mask(), "kmer is not masked"); }
    sized_kmer<kmer_t, ksize_t> make_sized(no_conversion_kmer<kmer_t> kmer) const
    {
        assert_masked(kmer);
        return {kmer, *this};
    }

  private:
    ksize_t _k;
};

/// Same as kmer_model but cache the mask
template<typename KmerT = size_t, typename KSizeT = ksize_t>
struct packed_layout kmer_mask : public kmer_model<KmerT, KSizeT>
{
    using base    = kmer_model<KmerT, KSizeT>;
    using kmer_t  = KmerT;
    using ksize_t = KSizeT;

    kmer_mask(kmer_model<kmer_t, ksize_t> k)
      : base(k)
      , _mask(base::mask())
    {}

    CPP14_CONSTEXPR kmer_t cardinality() const
    {
        assume(this->bitwidth() < bits::bitwidth<kmer_t>(), "cadinality is not representable in kmer type");
        return kmer_t(1) << this->bitwidth();
    }
    CPP14_CONSTEXPR kmer_t mask() const { return _mask; }
    CPP14_CONSTEXPR kmer_t mask(no_conversion_kmer<kmer_t> kmer) const { return kmer & mask(); }
    // void assume_masked(no_conversion_kmer<kmer_t> kmer) const { assume(kmer <= mask(), "kmer is not masked"); }
    void assert_masked(no_conversion_kmer<kmer_t> kmer) const { assert(kmer <= mask(), "kmer is not masked"); }
    sized_kmer<kmer_t, ksize_t> make_sized(no_conversion_kmer<kmer_t> kmer) const
    {
        assert_masked(kmer);
        return {kmer, *this};
    }

  private:
    kmer_t _mask;
};

template<typename kmer_t, typename ksize_t> struct packed_layout sized_kmer
{
    operator kmer_t() const { return kmer; }

    kmer_t                      kmer;
    kmer_model<kmer_t, ksize_t> size;
};

template<typename T, typename ksize_t = ksize_t>
constexpr inline auto
get_kmer(sized_kmer<T, ksize_t> sz_kmer) -> decltype(get_kmer(sz_kmer.kmer))
{
    return get_kmer(sz_kmer.kmer);
}

template<typename kmer_t = kmer_t, typename ksize_t = ksize_t>
std::string
to_string(const sized_kmer<kmer_t, ksize_t>& kmer)
{
    return kmer2str(kmer.kmer, kmer.size);
}

template<typename KmerT, typename From = KmerT, bool masking = true> struct sub_kextractor
{
    using kmer_t       = KmerT;
    using sized_kmer_t = sized_kmer<kmer_t>;

    sub_kextractor(ksize_t k, ksize_t offset)
      : mask_(bits::bitmask<From>(2 * k, 2 * offset))
      , shift_(2 * offset)
      , k_(k)
    {
        assert(2 * (k + offset) <= CHAR_BIT * sizeof(From), "k too large");
    }

    ksize_t size() const { return k_; }

    ksize_t bits() const { return 2 * k_; }

    size_t image_size() const { return size_t(mask_ >> shift_) + 1; }

    sized_kmer_t operator()(no_conversion_kmer<From> kmer) const { k_.make_sized(kmer); }

    ssize_t compare(no_conversion_kmer<From> a, no_conversion_kmer<From> b) const { return ssize_t(mask(a)) - mask(b); }

  protected:
    From mask(From x) const
    {
        if (masking) {
            return x & mask_;
        } else {
            assert(x <= this->mask_, "kmer larger than max value");
            return x;
        }
    }

    From               mask_;
    ksize_t            shift_;
    kmer_model<kmer_t> k_;
};

template<typename KmerT, typename From, bool masking = true> struct suffix_kextractor
{
    using kmer_t       = KmerT;
    using sized_kmer_t = sized_kmer<kmer_t>;

    suffix_kextractor(ksize_t k)
      : mask_(bits::bitmask<kmer_t>(2 * k))
      , k_(k)
    {
        assert(2 * k <= CHAR_BIT * sizeof(kmer_t), "k too large");
    }

    ksize_t size() const { return k_; }

    ksize_t bits() const { return k_.bits(); }

    size_t image_size() const { return size_t(mask_) + 1; }

    sized_kmer_t operator()(no_conversion_kmer<From> kmer) const { return {mask(kmer), k_}; }

    ssize_t compare(no_conversion_kmer<From> a, no_conversion_kmer<From> b) const
    {
        return typename std::make_signed<From>::type(mask(a)) - mask(b);
    }

  protected:
    kmer_t mask(From x) const
    {
        if (masking) {
            return x & mask_;
        } else {
            assert(x <= this->mask_, "kmer larger than max value");
            return x;
        }
    }

    kmer_t             mask_;
    kmer_model<kmer_t> k_;
};

template<typename KmerT, typename From, bool masking = false> struct prefix_kextractor
{
    using kmer_t = KmerT;

    prefix_kextractor(ksize_t k, ksize_t offset)
      : mask_(bits::bitmask<kmer_t>(2 * k))
      , shift_(2 * offset)
      , k_(k)
    {
        assert(2 * (k + offset) <= CHAR_BIT * sizeof(From), "k too large");
    }

    ksize_t size() const { return k_; }

    ksize_t bits() const { return k_.bits(); }

    size_t image_size() const { return mask_ + 1; }

    sized_kmer<kmer_t> operator()(no_conversion_kmer<From> kmer) const { return k_.make_sized(mask(kmer >> shift_)); }

    ssize_t compare(no_conversion_kmer<From> a, no_conversion_kmer<From> b) const
    {
        // FIXME: performance
        return ssize_t(mask(a >> shift_)) - ssize_t(mask(b >> shift_));
    }

  protected:
    kmer_t mask(From x) const
    {
        auto y = kmer_t(x);
        if (masking) {
            return y & mask_;
        } else {
            k_.assert_masked(y);
            assert(x <= this->mask_, "kmer larger than max value");
            return x;
        }
    }

    KmerT              mask_;
    ksize_t            shift_;
    kmer_model<kmer_t> k_;
};

namespace details {

/// Allows to compare things under a mapping object with a compare(x,y) function
template<typename Extractor, typename T> struct extracted_comparator
{
    const Extractor& ex;
    bool             operator()(const T& a, const T& b) { return this->ex.compare(a, b) < 0; }
};
} // namespace details

/// k-mer with a position in sequence
/// Defines an ordiering on the k-mer first, then the position
template<typename T> struct packed_layout positioned
{
    T      data;
    size_t pos;

    operator T() const { return data; }

    bool operator<(const positioned& other) const
    {
        ssize_t delta = ssize_t(data) - ssize_t(other.data);
        return likely(delta != 0) ? delta < 0 : this->pos < other.pos;
    }

    template<typename Extractor>
    static details::extracted_comparator<Extractor, positioned<T>> get_comparator(const Extractor& ex)
    {
        return {ex};
    }
};

template<typename T>
constexpr inline auto
get_kmer(positioned<T> pos_kmer) -> decltype(get_kmer(pos_kmer.data))
{
    return get_kmer(pos_kmer.data);
}

namespace details {
template<typename Extractor, typename T> struct extracted_comparator<Extractor, positioned<T>>
{
    const Extractor& ex;
    bool             operator()(const positioned<T>& a, const positioned<T>& b)
    {
        auto delta = ex.compare(a.data, b.data);
        static_assert(std::is_signed<decltype(delta)>::value, "unsigned substraction result");
        return likely(delta != 0) ? delta < 0 : a.pos < b.pos;
    }
};
} // namespace details

namespace details {

// TODO refine these helper as view facade for reuse in iterator_pair and friends
// Completes the definition of random iterator wrappers
template<typename Base> struct random_iter_wrapper : public Base
{
    using difference_type   = typename Base::difference_type;
    using reference         = typename Base::reference;
    using iterator_category = std::random_access_iterator_tag;
    using Base::Base;

    random_iter_wrapper(const random_iter_wrapper&) = default;
    random_iter_wrapper& operator=(const random_iter_wrapper&) = default;
    random_iter_wrapper(random_iter_wrapper&&)                 = default;
    random_iter_wrapper& operator=(random_iter_wrapper&&) = default;

    reference operator*() const { return {_repr}; }
    reference operator[](difference_type d) const { return {_repr + d}; }
    //    random_iter_wrapper& operator+=(difference_type d) { _repr += d; return *this; }
    random_iter_wrapper  operator+(difference_type d) { return {_repr + d}; }
    random_iter_wrapper& operator++()
    {
        ++_repr;
        return *this;
    }
    random_iter_wrapper operator++(int)
    {
        auto cpy = *this;
        ++*this;
        return cpy;
    }
    random_iter_wrapper& operator-=(difference_type d)
    {
        _repr -= d;
        return *this;
    }
    random_iter_wrapper  operator-(difference_type d) { return {_repr - d}; }
    random_iter_wrapper& operator--()
    {
        _repr--;
        return *this;
    }
    random_iter_wrapper operator--(int)
    {
        auto cpy = *this;
        --*this;
        return cpy;
    }
    difference_type operator-(const random_iter_wrapper& other) const { return this->_repr - other._repr; }
    bool            operator<(const random_iter_wrapper& other) const { return this->_repr < other._repr; }
    bool            operator<=(const random_iter_wrapper& other) const { return this->_repr <= other._repr; }
    bool            operator>(const random_iter_wrapper& other) const { return this->_repr > other._repr; }
    bool            operator>=(const random_iter_wrapper& other) const { return this->_repr >= other._repr; }
    bool            operator!=(const random_iter_wrapper& other) const { return this->_repr != other._repr; }
    bool            operator==(const random_iter_wrapper& other) const { return this->_repr == other._repr; }

  protected:
    using Base::_repr;
};

struct const_dna_bitstring_iter_base
{
  public:
    using difference_type = ptrdiff_t;
    using value_type      = nuc_t;

    struct reference
    {
        operator nuc_t()
        {
            ksize_t offset = _ptr & 0b11u;
            auto    ptr    = reinterpret_cast<const uint8_t * restrict>(_ptr >> 2);
            return as_nuc((*ptr >> offset) & 0b11u);
        }

      protected:
        template<typename C> friend struct random_iter_wrapper;
        reference(uintptr_t ptr)
          : _ptr(ptr)
        {}
        uintptr_t _ptr;
    };

    template<typename It>
    const_dna_bitstring_iter_base(It it, ksize_t off = 0)
      : _repr((reinterpret_cast<uintptr_t>(&*(it)) << 2) + off)
    {}

  protected:
    const_dna_bitstring_iter_base(uintptr_t repr)
      : _repr(repr)
    {}
    uintptr_t _repr;
};

struct dna_bitstring_iter_base : public const_dna_bitstring_iter_base
{
    struct reference : public const_dna_bitstring_iter_base::reference
    {
        reference& operator=(nuc_t n)
        {
            ksize_t offset = _ptr & 0b11;
            auto    ptr    = reinterpret_cast<uint8_t * restrict>(_ptr >> 2);
            *ptr |= uint8_t(n) << offset;
            return *this;
        }

      protected:
        using const_dna_bitstring_iter_base::reference::reference;
    };

    dna_bitstring_iter_base(const const_dna_bitstring_iter_base& x)
      : const_dna_bitstring_iter_base(x)
    {}
    template<typename It>
    dna_bitstring_iter_base(It it, ksize_t off = 0)
      : const_dna_bitstring_iter_base(it, off)
    {}
};

template<typename _repr_t> struct dna_ascii_iter_base
{
    using difference_type = ptrdiff_t;
    using value_type      = nuc_t;
    struct reference
    {
        operator nuc_t() const { return nuc2int(*_ptr); }

        const reference& operator=(nuc_t n) const { *const_cast<char*>(this->_ptr) = "ACTG"[as_integer(n)]; }

      protected:
        template<typename C> friend struct random_iter_wrapper;
        reference(_repr_t ptr)
          : _ptr(ptr)
        {}
        const _repr_t _ptr;
    };

    dna_ascii_iter_base(_repr_t ptr)
      : _repr(ptr)
    {}

  protected:
    _repr_t _repr;
};

template<typename R, typename _iterator_wrapper> struct wrapped_range
{
  private:
    using range_t = remove_reference_t<R>;

  public:
    using iterator   = _iterator_wrapper;
    using value_type = typename iterator::value_type;
    using reference  = typename iterator::reference;
    using size_type  = size_t;

    iterator begin() const { return _iterator_wrapper(gatbl::begin(repr)); }
    iterator end() const { return _iterator_wrapper(gatbl::end(repr)); }
    iterator begin() { return _iterator_wrapper(gatbl::begin(repr)); }
    iterator end() { return _iterator_wrapper(gatbl::end(repr)); }

    reference operator[](size_type i) const
    {
        assume(repr.begin() + i < repr.end(), "Past the end index");
        return *const_iterator(repr.begin() + i);
    }

    size_type size() const { return end() - begin(); }

    template<typename _R>
    wrapped_range(_R&& r)
      : repr(std::forward<_R>(r))
    {}
    wrapped_range(const wrapped_range&) = default;
    wrapped_range(wrapped_range&&)      = default;
    wrapped_range& operator=(const wrapped_range&) = default;
    wrapped_range& operator=(wrapped_range&&) = default;

    operator range_t&() { return repr; }
    operator const range_t&() const { return repr; }

    R repr; // left public for brace initialization
};

} // namespace details

using dna_bitstring_iter                       = details::random_iter_wrapper<details::dna_bitstring_iter_base>;
template<typename It> using dna_ascii_iter     = details::random_iter_wrapper<details::dna_ascii_iter_base<It>>;
template<typename R> using dna_bitstring_range = details::wrapped_range<R, dna_bitstring_iter>;
template<typename R> using dna_ascii_range     = details::wrapped_range<R, dna_ascii_iter<iterator_t<R>>>;

/// Represent an initiated iteration of a window over a range
template<typename R, typename W> struct window_range
{
    using window_t         = remove_reference_t<W>;
    using inner_iterator_t = gatbl::iterator_t<R>;
    using sentinel         = gatbl::sentinel_t<R>;
    using value_type       = remove_reference_t<decltype(*std::declval<window_t>())>;

    template<typename _W, typename... Args>
    window_range(_W&& w, Args&&... args) // That whole story is kinda sad: I'm not decided if reference or inlines are
                                         // better in most case and
      : _inner_range(std::forward<Args>(args)...)
      , _win(std::forward<_W>(w))
      , _it(_win.fill(_inner_range.begin()))
      , _end(_inner_range.end())
    {}

    // FIXME: view facade needed: this again is a meatless wrapper with all the logic in the iterator:
    using const_iterator = class iterator
    {
        window_range& _r;
        iterator(window_range& r)
          : _r(r){};
        friend struct window_range;

      public:
        flatten_fun hot_fun forceinline_fun iterator& operator++()
        {
            auto& it = _r._it;
            ++it;
            if (it != _r._end) _r._win.push_back(*it);
            return *this;
        }
        bool operator!=(sentinel sentinel) const flatten_fun hot_fun
        {
            // assume(&other._r == &_r, "Compare iterators from different window ranges");
            return _r._it != sentinel;
        }
        bool operator!=(const iterator& other) const flatten_fun hot_fun
        {
            assume(&other._r == &_r, "Compare iterators from different window ranges");
            return _r._it != other._r._end;
        }
        auto            operator*() -> decltype(*(_r._win)) const { return *(_r._win); }
        window_t&       window() { return _r._win; }
        const window_t& window() const { return _r._win; }
    };

    iterator begin() { return {*this}; }
    //    sentinel_t end() { return _end; }
    iterator end() { return {*this}; }

    window_t&       window() { return _win; }
    const window_t& window() const { return _win; }

    inner_iterator_t&       inner_iterator() { return _it; }
    const inner_iterator_t& inner_iterator() const { return _it; }

  private:
    R                _inner_range; // Let the parameterization choose if this a reference or not
    W                _win;         // Same here
    inner_iterator_t _it;
    const sentinel   _end;
};

template<typename R, typename W>
window_range<R, W>
operator|(R&& range, W&& window)
{
    return window_range<R, W>(std::forward<W>(window), std::forward<R>(range));
}

struct LexicoCanonical
{
    template<typename T> static pure_fun T canonize_bidir(T x, T y) { return x < y ? x : y; }

    template<typename T> static pure_fun T canonize(T x, ksize_t k) { return canonize_bidir(x, rcb(x, k)); }

    template<typename T> static pure_fun T uncanonize(T x) { return x; }
};

inline pure_fun int
popcount(unsigned x)
{
    return __builtin_popcount(x);
}
inline pure_fun int
popcount(unsigned long x)
{
    return __builtin_popcountl(x);
}
inline pure_fun int
popcount(unsigned long long x)
{
    return __builtin_popcountll(x);
}

struct ParityCanonical
{
    template<typename T> static pure_fun T canonize_bidir(T x, T y) { return (popcount(x) & 1 ? x : y) >> 1; }

    template<typename T> static pure_fun T canonize(T x, ksize_t k) { return (popcount(x) & 1 ? x : rcb(x, k)) >> 1; }

    template<typename T> static pure_fun T uncanonize(T x) { return T(popcount(x) & 1 ? 0 : 1) | (x << 1); }
};

template<typename kmer_t = kmer_t, typename Canonical = LexicoCanonical> struct kmer_window : private Canonical
{
    using value_type = kmer_t;

    kmer_window(ksize_t k)
      : _mask(~kmer_t(0) >> (CHAR_BIT * sizeof(kmer_t) - 2 * k))
      , _k(k)
      , _left_bitpos(2 * (k - 1))
    {
        assume(k <= sizeof(kmer_t) * 4, "k too large");
    }

    void clear()
    {
        _reverse = 0;
        _forward = 0;
    }

    template<typename It> hot_fun It fill(It it)
    {
        unchecked_push_back(*it);
        for (ksize_t i = 1; i < _k; ++i) {
            unchecked_push_back(*++it);
        }
        check();
        return it;
    }

    kmer_window& set_forward(kmer_t kmer)
    {
        _forward = kmer;
        _reverse = rcb(kmer, _k);
        return check();
    }

    kmer_window& set_reverse(kmer_t kmer)
    {
        _reverse = kmer;
        _forward = rcb(kmer, _k);
        return check();
    }

    hot_fun kmer_window& push_back(nuc_t nuc)
    {
        unchecked_push_back(nuc);
        mask(true);
        return check();
    }

    hot_fun kmer_window& push_front(nuc_t nuc)
    {
        unchecked_push_front(nuc);
        mask(false);
        return check();
    }

    void unchecked_push_back(nuc_t nuc)
    {
        _forward = (_forward << 2) | as_integer<kmer_t>(nuc);
        _reverse = (_reverse >> 2) | (kmer_t(0b10 ^ as_integer(nuc)) << _left_bitpos);
    }

    void unchecked_push_front(nuc_t nuc)
    {
        _forward = (_forward >> 2) | (as_integer<kmer_t>(nuc) << _left_bitpos);
        _reverse = (_reverse << 2) | (0b10 ^ as_integer<kmer_t>(nuc));
    }

    void mask(bool forward = true)
    {
        if (forward)
            _forward &= _mask;
        else
            _reverse &= _mask;
    }

    kmer_window& check()
    {
        assume(_forward <= _mask, "forward kmer is greater than max value");
        assume(_reverse <= _mask, "reverse kmer is greater than max value");
        assert(_forward == rcb(_reverse, _k), "Reversed sequence don't match the forward sequence");
        return *this;
    }

    using sized_kmer_t = sized_kmer<kmer_t, ksize_t>;
    ksize_t      size() const { return _k; }
    sized_kmer_t forward() const { return {_forward, _k}; }
    sized_kmer_t reverse() const { return {_reverse, _k}; }
    sized_kmer_t canon() const { return {Canonical::canonize_bidir(_forward, _reverse), _k}; }

    sized_kmer_t operator*() { return canon(); }

  private:
    const kmer_t    _mask;
    kmer_t          _forward = {}, _reverse = {};
    const ksize_t   _k;
    bits::bitsize_t _left_bitpos;
};

template<typename Base, typename Functor>
struct map_window
  : public Base
  , private Functor
{
    //    using Base::Base;

    template<typename... Args, typename F>
    map_window(F&& functor, Args&&... args)
      : Base(std::forward<Args>(args)...)
      , Functor(std::forward<F>(functor))
    {}

    using value_type = decltype(std::declval<Functor>()(*std::declval<Base>()));

    value_type operator*() { return Functor::operator()(Base::operator*()); }
};

/// Compute 2-mer entropy inside kmer for complexity filtering
/// The log2 entropy ranges from 0 to 4
template<typename kmer_t = kmer_t> class entropy_filter
{

    using lktnum_t   = uint16_t; // scaled (fixed point) entropy of a signle 2-mer
    using integral_t = size_t;   // scaled entropy

    static constexpr double precision = std::numeric_limits<lktnum_t>::max() * 1.88416938536372; // * log(2)/exp(-1)
    std::unique_ptr<lktnum_t[]>
                     _xlogx_lkt; // Lookup table giving the entropy contribution of a signle 2-mer given it's number of occurrences
    const integral_t _threshold;
    const ksize_t    _n;

    integral_t hot_fun _entropy(kmer_t kmer) const
    {
        uint8_t counts[16] = {0};

        for (int i = 0; i < _n; i++) {
            counts[kmer & kmer_t(15u)]++;
            kmer >>= 2;
        }
        assume(kmer < 4, "kmer larger than expected"); // A single base should remain

        integral_t ent = 0;
        for (int i = 0; i < 16; i++) {
            assume(counts[i] <= _n, "count out of range");
            ent += _xlogx_lkt[counts[i]];
        }

        return ent;
    }

  public:
    cold_fun entropy_filter(ksize_t k, double threshold)
      : _threshold(threshold * precision)
      , _n(k - 1) // Number of 2-mers
    {
        // Tabulate the -p*log2(p) values
        _xlogx_lkt     = make_unique<lktnum_t[]>(_n + 1);
        _xlogx_lkt[0]  = 0;
        _xlogx_lkt[_n] = 0;
        for (int i = 1; i < _n; i++) {
            double p = i / double(_n);
            // Full range value encountered for -p*log2(p) = 0 => p_max = exp(-1)
            // E_max = -log(p_max) / log(2) * p_max = exp(-1) / log(2)
            // precision = lktnum / E_max =  lktnum * 1.88417
            _xlogx_lkt[i] = lktnum_t(-log2(p) * p * precision);
            assume(size_t(-log2(p) * p * precision) <= size_t(std::numeric_limits<lktnum_t>::max()),
                   "out of range lookup entry");
        }
    }

    entropy_filter(entropy_filter&&) noexcept = default;
    entropy_filter(const entropy_filter& from)
      : _threshold(from._threshold)
      , _n(from._n)
    {
        _xlogx_lkt = make_unique<lktnum_t[]>(_n + 1);
        std::copy(from._xlogx_lkt.get(), from._xlogx_lkt.get() + _n + 1, _xlogx_lkt.get());
    }

    double entropy(kmer_t kmer) const { return double(_entropy(kmer)) / precision; }

    integral_t numerator(kmer_t kmer) const { return _entropy(kmer); }

    static constexpr integral_t denominator(kmer_t = 0) { return precision; }

    bool hot_fun operator()(kmer_t kmer) const { return _entropy(kmer) >= _threshold; }
};

/// A functor returning a kmer along with it's hash with an ordering based on the hash
template<typename HashFunctor = ReversibleHash> struct hash_kmer : private HashFunctor
{
    using HashFunctor::HashFunctor;

    template<typename _kmer_t, typename hash_t> struct hashed_kmer
    {
        using kmer_t = remove_reference_t<_kmer_t>;
        kmer_t kmer;
        hash_t hash;
               operator kmer_t() const { return kmer; }
        bool   operator<=(const hashed_kmer& other) const { return hash <= other.hash; }
    };

    template<typename kmer_t> auto operator()(kmer_t&& x) -> hashed_kmer<kmer_t, decltype(HashFunctor::operator()(x))>
    {
        auto hash = HashFunctor::operator()(x);
        return {std::forward<kmer_t>(x), std::move(hash)};
    }
};

template<typename rank_t> struct rank_kmer
{
    rank_kmer(const rank_t* ranks)
      : _ranks(ranks){};

    template<typename _kmer_t> struct ranked_kmer
    {
        using kmer_t = remove_reference_t<_kmer_t>;
        kmer_t kmer;
        rank_t rank;
               operator kmer_t() const { return kmer; }
        bool   operator<=(const ranked_kmer& other) const { return rank <= other.rank; }
    };

    template<typename kmer_t> auto operator()(kmer_t x) -> ranked_kmer<kmer_t>
    {
        auto rank = _ranks[x];
        return {std::forward<kmer_t>(x), std::move(rank)};
    }

  protected:
    const rank_t* _ranks;
};

template<typename kmer_t = uint64_t, typename Canonical = LexicoCanonical>
window_range<dna_ascii_range<const std::string&>, kmer_window<kmer_t, Canonical>>
iter_kmers(const std::string& str, ksize_t k)
{
    return {str, k};
};

template<typename kmer_t = uint64_t, typename Canonical = LexicoCanonical, typename R>
window_range<R, kmer_window<kmer_t, Canonical>>
iter_kmers(R str, ksize_t k)
{
    return {std::forward<R>(str), k};
};

template<typename elem_t, size_t log2size = 5, typename idx_t = ksize_t> struct stack_ringbuf
{
    static constexpr size_t capacity = 1ull << log2size;
    using value_type                 = elem_t;
    static_assert(capacity <= std::numeric_limits<idx_t>::max(), "idx_t too short");

  private:
    static constexpr size_t mask = capacity - 1;
    struct state
    {
        idx_t _b = 0;
        idx_t _e = 0;
    };

    template<typename ref_t = const stack_ringbuf&> struct accessor_tmpl : protected state
    {
        const value_type& front() const
        {
            check_not_empty();
            return at(state::_b);
        }
        const value_type& back() const
        {
            check_not_empty();
            return at(state::_e - 1);
        }

        bool  empty() const { return state::_b == state::_e; }
        bool  full() const { return size() == capacity; }
        idx_t size() const
        {
            assume(idx_t(state::_e - state::_b) < capacity, "Unwrapped indices");
            return state::_e - state::_b;
        }

      protected:
        friend class stack_ringbuf;
        accessor_tmpl(ref_t ref)
          : _ref(ref)
          , state(ref._state)
        {}

        void check_not_full() const { assume(!full(), "Ring full"); }
        void check_not_empty() const { assume(!empty(), "Ring empty"); }

        const value_type& at(idx_t idx) const
        {
            idx &= mask;
            assume(idx < capacity, "WTF");
            return _ref._arr[idx & mask];
        }
        ref_t _ref;
    };

    struct mutator : public accessor_tmpl<stack_ringbuf&>
    {
        // Here's the trick: copy back the cursors once we're done with mutations:
        ~mutator() { this->_ref._state = *this; }

        void clear()
        {
            // state::_e = state::_b;
            state::_e = state::_b = 0;
        }

        value_type& push_front(value_type val)
        {
            this->check_not_full();
            return this->at(--state::_b) = val;
        }
        value_type& push_back(value_type val)
        {
            this->check_not_full();
            return this->at(state::_e++) = val;
        }
        value_type& front()
        {
            this->check_not_empty();
            return this->at(state::_b);
        }
        value_type& back()
        {
            this->check_not_empty();
            return this->at(state::_e - 1);
        }
        value_type pop_front()
        {
            this->check_not_empty();
            return std::move(this->at(state::_b++));
        }
        value_type pop_back()
        {
            this->check_not_empty();
            return std::move(this->at(--state::_e));
        }

      protected:
        using accessor_tmpl<stack_ringbuf&>::accessor_tmpl;
        value_type& at(idx_t idx)
        {
            idx &= mask;
            assume(idx < capacity, "WTF");
            return this->_ref._arr[idx & mask];
        }
    };

    elem_t _arr[capacity];
    state  _state;

  public:
    using accessor = const accessor_tmpl<const stack_ringbuf&>;

    accessor operator*() const { return {*this}; }
    mutator  operator*() { return {*this}; }
};

template<typename elem_t, size_t log2size = 5, typename idx_t = ksize_t> struct stack_ringbuf_old
{
    static constexpr size_t capacity = 1ull << log2size;
    static_assert(capacity <= std::numeric_limits<idx_t>::max(), "idx_t too short");
    elem_t& push_front(elem_t val)
    {
        check_not_full();
        return _arr[mask & --_b] = val;
    }
    elem_t& push_back(elem_t val)
    {
        check_not_full();
        return _arr[mask & _e++] = val;
    }
    elem_t& front()
    {
        check_not_empty();
        return _arr[mask & _b];
    }
    elem_t& back()
    {
        check_not_empty();
        return _arr[mask & (_e - 1)];
    }
    const elem_t& front() const
    {
        check_not_empty();
        return _arr[mask & _b];
    }
    const elem_t& back() const
    {
        check_not_empty();
        return _arr[mask & (_e - 1)];
    }
    elem_t pop_front()
    {
        check_not_empty();
        return std::move(_arr[mask & _b++]);
    }
    elem_t pop_back()
    {
        check_not_empty();
        return std::move(_arr[mask & --_e]);
    }
    void  clear() { _e = _b; }
    bool  empty() const { return _b == _e; }
    bool  full() const { return size() == capacity; }
    idx_t size() const
    {
        assume(_e - _b < capacity, "Unwrapped indices");
        return _e - _b;
    }

    const stack_ringbuf_old& operator*() const { return {*this}; }
    stack_ringbuf_old&       operator*() { return {*this}; }

  private:
    void check_not_full() const { assume(!full(), "Ring full"); }
    void check_not_empty() const { assume(!empty(), "Ring empty"); }

    static constexpr size_t mask = capacity - 1;
    elem_t                  _arr[capacity];
    idx_t                   _b = 0;
    idx_t                   _e = 0;
};

template<typename elem_t, size_t log2size = 5, typename idx_t = ksize_t> struct heap_ringbuf
{
    static constexpr size_t capacity = 1ull << log2size;
    static_assert(capacity <= std::numeric_limits<idx_t>::max(), "idx_t too short");

    heap_ringbuf()
      : _arr(new elem_t[capacity])
    {}
    heap_ringbuf(const heap_ringbuf&) = delete;
    heap_ringbuf(heap_ringbuf&&)      = default;

    elem_t& push_front(elem_t val)
    {
        check_not_full();
        return _arr[mask & --_b] = val;
    }
    elem_t& push_back(elem_t val)
    {
        check_not_full();
        return _arr[mask & _e++] = val;
    }
    elem_t& front()
    {
        check_not_empty();
        return _arr[mask & _b];
    }
    elem_t& back()
    {
        check_not_empty();
        return _arr[mask & (_e - 1)];
    }
    const elem_t& front() const
    {
        check_not_empty();
        return _arr[mask & _b];
    }
    const elem_t& back() const
    {
        check_not_empty();
        return _arr[mask & (_e - 1)];
    }
    elem_t pop_front()
    {
        check_not_empty();
        return std::move(_arr[mask & _b++]);
    }
    elem_t pop_back()
    {
        check_not_empty();
        return std::move(_arr[mask & --_e]);
    }
    void  clear() { _e = _b; }
    bool  empty() const { return _b == _e; }
    bool  full() const { return size() == capacity; }
    idx_t size() const
    {
        assume(_e - _b < capacity, "Unwrapped indices");
        return _e - _b;
    }

    const heap_ringbuf& operator*() const { return {*this}; }
    heap_ringbuf&       operator*() { return {*this}; }

  private:
    void check_not_full() const { assume(!full(), "Ring full"); }
    void check_not_empty() const { assume(!empty(), "Ring empty"); }

    static constexpr size_t   mask = capacity - 1;
    std::unique_ptr<elem_t[]> _arr;
    idx_t                     _b = 0;
    idx_t                     _e = 0;
};

template<typename elem_t = minimizer_t, bits::bitsize_t log2maxsize = 5, typename idx_t = size_t> struct minimum_window
{
    minimum_window(ksize_t w)
      : _w(w)
    {}

    idx_t size() const { return _w; }

    struct value_type
    {
        unsigned pos;
        elem_t   elem;
                 operator elem_t&() { return elem; }
                 operator const elem_t&() const { return elem; }
    };

    static_assert(std::is_pod<elem_t>::value, "not a POD");

    template<typename It> It fill(It it)
    {
        auto& ring = *_ring;
        ring.clear();

        value_type rec;
        rec.pos  = 0;
        rec.elem = *it;
        ring.push_back(rec);

        // Fill the first window
        for (; rec.pos < _w; rec.pos++) {
            elem_t elem = *++it;

            if (elem <= rec.elem) {
                do
                    ring.pop_back();
                while (not ring.empty() && elem <= ring.back().elem);
            }

            rec.elem = elem;
            ring.push_back(rec);
        }
        return it;
    }

    hot_fun bool push_back(elem_t elem)
    {
        auto&      ring = *_ring;
        value_type rec  = ring.back();

        bool empty = false;
        if (elem <= rec.elem) {
            do
                ring.pop_back();
            while (not ring.empty() && elem <= ring.back().elem);
            empty = ring.empty();
        }

        rec.pos++;
        rec.elem = elem;
        ring.push_back(rec);

        if (not empty) { // No new minimum entering the right side
            // Check for the current minimum going out on the left
            if (ring.front().pos + _w > rec.pos) {
                return false;
            } else {
                ring.pop_front();
                return true;
            }
        } else {
            return true;
        }
    }

    const value_type& operator*() const { return (*_ring).front(); }

    stack_ringbuf_old<value_type, log2maxsize, idx_t> _ring;
    const idx_t                                       _w;
};

template<typename elem_t = minimizer_t, bits::bitsize_t log2maxsize = 5> struct minimum_window2
{
    minimum_window2(ksize_t w)
      : _w(w)
    {}

    ksize_t size() const { return _w; }

    struct value_type
    {
        unsigned pos;
        elem_t   elem;
                 operator elem_t&() { return elem; }
                 operator const elem_t&() const { return elem; }
    };

    static_assert(std::is_pod<elem_t>::value, "not a POD");

    template<typename It> It fill(It it)
    {
        _ring.clear();

        value_type rec;
        rec.pos  = 0;
        rec.elem = *it;
        front    = rec;
        //_ring.push_back(rec);

        // Fill the first window
        for (; rec.pos < _w; rec.pos++) {
            elem_t elem = *++it;
            rec.elem    = *++it;

            if (rec.elem <= front) {
                front = rec;
                _ring.clear();
            } else {
                while (not _ring.empty() && rec.elem <= _ring.back().elem)
                    _ring.pop_back();
                if (not _ring.empty()) {
                    _ring.push_back(rec);
                } else {
                    front = rec;
                }
            }
        }

        return it;
    }

    bool push_back(elem_t elem)
    {
        value_type rec;
        if (not _ring.empty()) {
            value_type rec = _ring.back();

            // 93593784
            bool empty = false;
            if (elem <= rec.elem) {
                // 46925816 = 50%
                do
                    _ring.pop_back();
                while (not _ring.empty() && elem <= _ring.back().elem);
                empty = _ring.empty();
            }
        } else {
            rec = front;
        }

        rec.pos++;
        rec.elem = elem;
        _ring.push_back(rec);

        if (not _ring.empty()) { // No new minimum entering the right side
            // Check for the current minimum going out on the left
            if (_ring.front().pos + _w > rec.pos) {
                return false;
            } else {
                _ring.pop_front();
                return true;
            }
            // assume(_ring.front().pos > rec.pos - _w, "WTF");
        } else {
            return true;
        }
    }

    const value_type& operator*() const { return _ring.front(); }

    value_type                             front;
    stack_ringbuf<value_type, log2maxsize> _ring;
    ksize_t                                _w;
};

} // namesapce gatbl
#endif // KMER_HPP
