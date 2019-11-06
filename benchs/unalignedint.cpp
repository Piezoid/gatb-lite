#include <vector>

#include <benchmark/benchmark.h>
#include "gatbl/ds/int_vector.hpp"

template<typename T>
inline void
use(T&& t)
{
    __asm__ __volatile__("" ::"g"(t));
}

void
clobber()
{
    __asm__ __volatile__("" : : : "memory");
}

inline uint64_t
splitmix64_stateless(uint64_t index)
{
    uint64_t z = (index + UINT64_C(0x9E3779B97F4A7C15));
    z          = (z ^ (z >> 30u)) * UINT64_C(0xBF58476D1CE4E5B9);
    z          = (z ^ (z >> 27u)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31u);
}

template<typename T = uint64_t>
static inline std::vector<T>
make_randint_vector(size_t n, uint64_t bits, size_t seed = 1)
{
    std::vector<T> res(n);
    const T        mask = T((size_t(1) << bits) - 1);
    for (size_t i = 0; i < n; i++)
        res[i] = T(splitmix64_stateless(i * seed)) & mask;
    return res;
}

template<typename T = uint64_t>
static inline std::vector<T>
make_sorted_randint_vector(size_t n, uint64_t bits, size_t seed = 1)
{
    auto res = make_randint_vector<T>(n, bits, seed);
    std::sort(res.begin(), res.end());
    return res;
}

constexpr inline size_t
gcd(size_t __m, size_t __n)
{
    return __m == 0 ? __n : __n == 0 ? __m : gcd(__n, __m % __n);
}

constexpr inline size_t
lcm(size_t __m, size_t __n)
{
    return (__m != 0 && __n != 0) ? (__m / gcd(__m, __n)) * __n : 0;
}

using gatbl::bits::bitsize_t;

struct unaligned_ref
{
  public:
    using type = size_t;

  protected:
    using ptr_ty                     = uint8_t*;
    static constexpr size_t bitwidth = sizeof(type) * CHAR_BIT;

    type mask() const noexcept { return type(-1) >> (bitwidth - len); }

    unaligned_ref(ptr_ty p, bitsize_t o, bitsize_t w) noexcept
      : ptr(p)
      , len(w)
      , offset(o)
    {
        assume(len + offset <= bitwidth, "word overflow");
    }

    ptr_ty    ptr;
    bitsize_t len, offset;

  public:
    operator type() const noexcept
    {
        type w;
        memcpy(&w, ptr, sizeof(type));
        return (w >> offset) & mask();

        w <<= uint8_t(64u - len - offset);
        w >>= uint8_t(64u - len);
        return w;
    }
    type get() const noexcept { return *this; }

    const unaligned_ref& operator=(type x) const noexcept
    {
        type w;
        memcpy(&w, ptr, sizeof(type));

        type m = mask() << offset;
        assume(x <= mask(), "value out of bounds");
        w = (w & ~m) | (x << offset);

        memcpy(ptr, &w, sizeof(type));
        return *this;
    }

    static unaligned_ref make_from_pos(type* baseptr, size_t pos, bitsize_t width)
    {
        size_t    bitoffset = pos * width;
        auto      ptr       = reinterpret_cast<ptr_ty>(baseptr) + (bitoffset >> 3);
        bitsize_t offset    = bitsize_t(7) & bitoffset;
        return {ptr, offset, width};
    }

    static unaligned_ref make_from_bitoffset(type* baseptr, size_t bitoffset, bitsize_t width)
    {
        auto      ptr    = reinterpret_cast<ptr_ty>(baseptr) + (bitoffset >> 3);
        bitsize_t offset = bitsize_t(7) & bitoffset;
        return {ptr, offset, width};
    }

    static unaligned_ref make_from_bitaddr(ptrdiff_t bitaddr, bitsize_t width)
    {
        return {reinterpret_cast<uint8_t*>(bitaddr >> 3), bitsize_t(bitaddr & 7u), width};
    }

    static unaligned_ref make_from_ptr_and_offset(ptr_ty ptr, bitsize_t offset, bitsize_t width)
    {
        return {ptr, offset, width};
    }

    friend void swap(unaligned_ref a, unaligned_ref b) noexcept
    { // FIXME: masked XOR swap may be faster
        type tmp = a;
        a        = type(b);
        b        = tmp;
    }
};

size_t
get(const unaligned_ref& r) noexcept
{
    return r;
}

struct cachemask_unalignedref : public unaligned_ref
{
    using unaligned_ref::type;

    operator type() const noexcept
    {
        type w;
        memcpy(&w, ptr, sizeof(type));
        return (w >> offset);

        w <<= uint8_t(64u - len - offset);
        w >>= uint8_t(64u - len);
        return w;
    }

    const unaligned_ref& operator=(type x) const noexcept
    {
        type w;
        memcpy(&w, ptr, sizeof(type));

        type m = mask << offset;
        assume(x <= mask, "value out of bounds");
        w = (w & ~m) | (x << offset);

        memcpy(ptr, &w, sizeof(type));
        return *this;
    }

    cachemask_unalignedref(unaligned_ref&& ref)
      : unaligned_ref(std::move(ref))
      , mask(unaligned_ref::mask())
    {}

  protected:
    type mask;
};

struct stepped_masked_iterator
{
    using value_type = size_t;

    stepped_masked_iterator() noexcept = default;

    stepped_masked_iterator(value_type* p, bitsize_t w, bitsize_t s)
      : addr(p)
      , width(w)
      , step(s)
    {}

    stepped_masked_iterator operator+(ptrdiff_t x) const
    {
        stepped_masked_iterator r = *this;
        r += x;
        return r;
    }

    stepped_masked_iterator& operator+=(ptrdiff_t x)
    {
        addr += x * step;
        return *this;
    }

    stepped_masked_iterator& operator++()
    {
        addr += step;
        return *this;
    }

    ssize_t operator-(const stepped_masked_iterator& other) const
    {
        auto d = (addr - other.addr) / step;
        assume(d % step == 0, "unaligned iterator");
        return d / step;
    }

    bool operator==(const stepped_masked_iterator& other) const { return addr == other.addr; }
    bool operator!=(const stepped_masked_iterator& other) const { return addr != other.addr; }
    bool operator<(const stepped_masked_iterator& other) const { return addr < other.addr; }

    struct reference
    {
        using type = stepped_masked_iterator::value_type;

        operator type() const noexcept { return *ptr & mask; }

        reference operator=(type x) const noexcept
        {
            assume(x <= mask, "int larger than expected");
            *ptr = (*ptr & ~mask) | x;
            return *this;
        }

      private:
        friend struct stepped_masked_iterator;
        static constexpr size_t bitwidth = sizeof(type) * CHAR_BIT;
        reference(value_type* p, value_type width) noexcept
          : ptr(p)
          , mask(value_type(-1) >> (bitwidth - width))
        {
            assume(width <= bitwidth, "word overflow");
        }
        value_type* ptr;
        value_type  mask;
    };

    reference operator*() const { return {addr, width}; }
    reference operator[](size_t i) const noexcept { return *(*this + ptrdiff_t(i)); }

  protected:
    value_type* addr; // A signed type enable arithmetic right shift for restoring a cannonical x86_64 pointer
    bitsize_t   width;
    bitsize_t   step;
};

struct it_ptrtrick
  : public gatbl::
      iterator_facade<it_ptrtrick, std::random_access_iterator_tag, size_t, unaligned_ref, void, std::ptrdiff_t>
{
  protected:
    using base = gatbl::
      iterator_facade<it_ptrtrick, std::random_access_iterator_tag, size_t, unaligned_ref, void, std::ptrdiff_t>;

  public:
    //    it_ptrtrick() noexcept = default;

    it_ptrtrick(value_type* p, bitsize_t w, size_t i = 0)
      : addr((reinterpret_cast<ptrdiff_t>(p) << 3) + ptrdiff_t(i * w))
      , width(w)
    {}

    it_ptrtrick& operator+=(ptrdiff_t x)
    {
        addr += x * width;
        return *this;
    }

    friend ssize_t compare(const it_ptrtrick& lhs, const it_ptrtrick& rhs) { return lhs.addr - rhs.addr; }

    using base::operator-;
    ssize_t     operator-(const it_ptrtrick& other) const
    {
        ptrdiff_t d = addr - other.addr;
        assume(d % width == 0, "unaligned addresses compared");
        return (addr - other.addr) / width;
    }

    unaligned_ref operator*() const { return reference::make_from_bitaddr(addr, width); }

  protected:
    ptrdiff_t addr; // A signed type enable arithmetic right shift for restoring a cannonical x86_64 pointer
    bitsize_t width;
};

template<typename T>
struct integral_iterator
  : gatbl::iterator_facade<integral_iterator<T>, std::random_access_iterator_tag, T, T, T*, gatbl::make_signed_t<T>>
{
    using difference_type = gatbl::make_signed_t<T>;
    T&                 operator*() noexcept { return _value; }
    T&                 operator*() const noexcept { return _value; }
    difference_type    operator-(T other) const noexcept { return _value - other; }
    difference_type    operator-(integral_iterator other) const noexcept { return _value - other._value; }
    integral_iterator& operator+=(difference_type d) noexcept
    {
        _value += d;
        return *this;
    }

    integral_iterator(T value)
      : _value(value)
    {}
    operator T() const noexcept { return _value; }

  protected:
    T _value;
};

struct it_bitoffset
  : gatbl::iterator_facade<it_bitoffset, std::random_access_iterator_tag, size_t, unaligned_ref, void, ssize_t>
{
  protected:
    using base
      = gatbl::iterator_facade<it_bitoffset, std::random_access_iterator_tag, size_t, unaligned_ref, void, ssize_t>;

  public:
    //    it_bitoffset() noexcept = default;

    it_bitoffset(value_type* p, bitsize_t w, size_t i = 0)
      : baseaddr(p)
      , bitoffset(i * w)
      , width(w)
    {}

    it_bitoffset& operator+=(ssize_t x)
    {
        ssize_t tmp = ssize_t(bitoffset) + (x * width);
        assume(tmp >= 0, "out of bound iterator");
        bitoffset = size_t(tmp);
        return *this;
    }

    friend ssize_t compare(const it_bitoffset& lhs, const it_bitoffset& rhs)
    {
        assume(lhs.baseaddr == rhs.baseaddr, "unconsistend base address");
        return ssize_t(lhs.bitoffset) - ssize_t(rhs.bitoffset);
    }

    using base::operator-;
    ssize_t     operator-(const it_bitoffset& other) const
    {
        ssize_t d = compare(*this, other);
        assume(d % width == 0, "unaligned address compared");
        return d / width;
    }

    reference operator*() const { return reference::make_from_bitoffset(baseaddr, bitoffset, width); }

  protected:
    value_type* baseaddr;
    size_t      bitoffset;
    bitsize_t   width;
};

struct it_pos : gatbl::iterator_facade<it_pos, std::random_access_iterator_tag, size_t, unaligned_ref, void, ssize_t>
{
  protected:
    using base = gatbl::iterator_facade<it_pos, std::random_access_iterator_tag, size_t, unaligned_ref, void, ssize_t>;

  public:
    //    it_pos() noexcept = default;

    it_pos(value_type* p, bitsize_t w, size_t i = 0)
      : baseaddr(p)
      , pos(i)
      , width(w)
    {}

    it_pos& operator+=(ssize_t x)
    {
        ssize_t new_pos = ssize_t(pos) + x;
        assume(new_pos >= 0, "out of range iterator");
        pos = size_t(new_pos);
        return *this;
    }

    using base::operator-;
    ssize_t     operator-(const it_pos& other) const
    {
        check_consistent_baseaddr(other);
        return ssize_t(pos) - ssize_t(other.pos);
    }

    reference operator*() const
    {
        return reference::make_from_bitoffset(const_cast<size_t* __restrict>(baseaddr), pos * width, width);
    }
    reference operator[](size_t i) const noexcept { return *(*this + ssize_t(i)); }

  protected:
    void check_consistent_baseaddr(const it_pos& other) const
    {
        assume(baseaddr == other.baseaddr, "unconsistend base address");
    }

    value_type* baseaddr;
    size_t      pos;
    bitsize_t   width;
};

template<typename T>
static inline T
refload(const T* __restrict _p, size_t i, bitsize_t n) noexcept
{
    auto* p = reinterpret_cast<const uint8_t*>(_p);
    if (n > 64 - 8) __builtin_unreachable();

    uint8_t mask_zeros = 64 - n;

    size_t    first_bit = i * n;
    bitsize_t idx_bits  = first_bit & 7u;
    size_t    idx_bytes = first_bit >> 3u;
    T         mask      = T(-1) >> (mask_zeros);

    T w;
    memcpy(&w, p + idx_bytes, sizeof(w));
    w >>= idx_bits;
    w &= mask;
    if (w >= (T(1) << n)) __builtin_unreachable();
    return w;
}

template<typename It = it_ptrtrick> struct int_mockup
{
    using iterator  = It;
    using reference = typename iterator::reference;
    using block_t   = typename iterator::value_type;

    static constexpr size_t block_size = sizeof(block_t);
    static constexpr size_t block_bits = block_size * CHAR_BIT;

    // NB: we loose block_bits-width on the first word since the first int is aligned to the low bits of the first word
    int_mockup(size_t sz, bitsize_t width)
      : memory((sz * width + (block_bits - width) + block_bits - 1) / block_bits)
      , _size(sz)
      , mask(block_t(-1) >> (block_bits - width))
      , width(width)
    {
        auto gcd   = ::gcd(block_bits, width);
        int_sync   = block_bits / gcd;
        block_sync = width / gcd;
        assert(size() >= sz, "oops");
    }

    size_t size() const noexcept { return _size; }

    template<typename T>
    int_mockup(std::vector<T> vec, T max_value)
      : int_mockup(vec.size(), bitsize_t(gatbl::bits::ilog2p1(max_value)))
    {
        const auto last = end();
        {
            auto it = begin();
            for (size_t i = 0; i < vec.size(); i++, ++it) {
                assume(it < last, "last reached");
                *it = block_t(vec[i]);
                assert(*it == block_t(vec[i]), "wtf");
            }
        }
#ifndef NDEBUG
        {
            auto it = begin();
            for (size_t i = 0; i < vec.size(); i++, ++it) {
                assume(it < last, "last reached");
                assert(*it == block_t(vec[i]), "wtf");
                assert(refload(memory.data(), i, width) == block_t(vec[i]), "wtf");
            }
        }

        {
            auto       it   = stepped_masked_iterator(memory.data(), width, block_sync);
            const auto last = stepped_masked_iterator(&*memory.end(), width, block_sync);
            for (size_t i = 0; i * int_sync < vec.size(); i += int_sync, ++it) {
                assume(it < last, "last reached");
                assert(*it == block_t(vec[i]), "wtf");
            }
        }
#endif
    }

    iterator begin() { return {memory.data(), width}; }
    iterator end() { return {memory.data(), width, _size}; }

    reference operator[](size_t i) noexcept { return *(begin() + ptrdiff_t(i)); }
    //    block_t operator[](size_t i) const noexcept { return *(const_cast<int_mockup*>(this)->begin() + i); }
    block_t operator[](size_t i) const noexcept { return *iterator(const_cast<block_t*>(memory.data()), width, i); }

    std::vector<block_t> memory;
    size_t               _size;
    block_t              mask;
    bitsize_t            width, block_sync, int_sync;
};

/// Baseline testing
struct int_mockup0 : int_mockup<>
{
    using int_mockup::int_mockup;

    block_t operator[](size_t i) const noexcept { return i; }
};

struct int_mockup1 : int_mockup<>
{
    using int_mockup::int_mockup;
    using typename int_mockup::block_t;

    block_t operator[](size_t i) const noexcept
    {
        block_t x = refload(memory.data(), i, width);
        assert(x == int_mockup::operator[](i), "wtf");
        return x;
    }
};

struct int_mockupref : int_mockup<>
{
    using int_mockup::int_mockup;
    using typename int_mockup::block_t;
    using typename int_mockup::reference;

    block_t operator[](size_t i) const noexcept
    {
        block_t x = load(memory.data(), i, width);
        assert(x == int_mockup::operator[](i), "wtf");
        return x;
    }

    template<typename T> static block_t load(const T* __restrict _p, size_t i, bitsize_t n) noexcept
    {
        return reference::make_from_bitoffset(const_cast<size_t* __restrict>(_p), i * n, n);
    }
};

template<typename Vec>
static noinline_fun void
do_query(const Vec& v, size_t mask, uint64_t seed, size_t nqueries) noexcept
{
    size_t x = 0;
    for (size_t i = 0; i < nqueries; i++) {
        size_t idx = splitmix64_stateless(i * (seed + 1));
        x += v[mask & idx]; // + v[mask & size_t(-idx)];
    }
    benchmark::DoNotOptimize(x);
}

template<typename Vec>
static void
random_access(benchmark::State& state)
{

    constexpr uint8_t  width       = 17;
    constexpr uint32_t max_element = (size_t(1) << width) - 1;
    constexpr size_t   log2        = 13;
    constexpr size_t   N           = size_t(1) << log2;
    constexpr size_t   idx_mask    = N - 1;

    constexpr size_t query_size = 1024;

    Vec v(make_sorted_randint_vector<uint32_t>(N, width), max_element);

    unsigned seed = 0;
    for (auto _ : state) {
        do_query(v, idx_mask, ++seed, query_size);
    }
}
BENCHMARK_TEMPLATE(random_access, int_mockup0);
BENCHMARK_TEMPLATE(random_access, int_mockup1);
BENCHMARK_TEMPLATE(random_access, int_mockupref);
BENCHMARK_TEMPLATE(random_access, int_mockup<it_ptrtrick>);
BENCHMARK_TEMPLATE(random_access, int_mockup<it_bitoffset>);
BENCHMARK_TEMPLATE(random_access, int_mockup<it_pos>);
BENCHMARK_TEMPLATE(random_access, gatbl::int_vector<>);

template<typename Vec>
static noinline_fun void
do_sort(Vec& v) noexcept
{
    std::sort(v.begin(), v.end());
    benchmark::DoNotOptimize(v);
}

static void
sort_stdvec(benchmark::State& state)
{

    constexpr uint8_t width = 7;
    constexpr size_t  log2  = 14;
    constexpr size_t  N     = size_t(1) << log2;

    auto vec = make_randint_vector<uint16_t>(N, width);

    for (auto _ : state) {
        state.PauseTiming();
        std::vector<uint16_t> v = vec;
        state.ResumeTiming();
        do_sort(v);
    }
}
BENCHMARK(sort_stdvec);

template<typename Vec>
static void
sort_intvec(benchmark::State& state)
{
    constexpr uint8_t  width       = 7;
    constexpr uint16_t max_element = (size_t(1) << width) - 1;
    constexpr size_t   log2        = 14;
    constexpr size_t   N           = size_t(1) << log2;

    auto vec = make_randint_vector<uint16_t>(N, width);

    for (auto _ : state) {
        state.PauseTiming();
        Vec v(vec, max_element);
        state.ResumeTiming();
        do_sort(v);
    }
}
BENCHMARK_TEMPLATE(sort_intvec, int_mockup<it_ptrtrick>);
BENCHMARK_TEMPLATE(sort_intvec, int_mockup<it_bitoffset>);
BENCHMARK_TEMPLATE(sort_intvec, int_mockup<it_pos>);
BENCHMARK_TEMPLATE(sort_intvec, gatbl::int_vector<>);

template<typename It = int* __restrict, typename T = gatbl::value_t<It>>
static inline void
binsearch4(It first, int howmany, size_t arraysize, T* __restrict targets, size_t* __restrict solutions)
{
    if (arraysize == 0) return; // pathological
    size_t i = 0;
    if (howmany >= 4) {
        for (; i < howmany; i += 4) {
            It     base0 = first;
            It     base1 = first;
            It     base2 = first;
            It     base3 = first;
            size_t n     = arraysize;
            while (n > 1) {
                size_t half = n >> 1;
                auto   mid0 = base0 + half;
                auto   mid1 = base1 + half;
                auto   mid2 = base2 + half;
                auto   mid3 = base3 + half;
                base0       = (*mid0 < targets[i + 0]) ? mid0 : base0;
                base1       = (*mid1 < targets[i + 1]) ? mid1 : base1;
                base2       = (*mid2 < targets[i + 2]) ? mid2 : base2;
                base3       = (*mid3 < targets[i + 3]) ? mid3 : base3;
                n -= half;
            }
            solutions[i + 0] = (*base0 < targets[i + 0]) + (base0 - first);
            solutions[i + 1] = (*base1 < targets[i + 1]) + (base1 - first);
            solutions[i + 2] = (*base2 < targets[i + 2]) + (base2 - first);
            solutions[i + 3] = (*base3 < targets[i + 3]) + (base3 - first);
        }
    }
    assume(howmany == i, "nope");
    //    for (; i < howmany; i++) {
    //        solutions[i] = branchfree_search(data[i], arraysize, targets[i]);
    //    }
}

template<int N = 4, typename Vec, typename It = gatbl::iterator_t<Vec>, typename T = gatbl::value_t<It>>
static inline void
binsearch(Vec& v, int howmany, T* __restrict targets, size_t* __restrict solutions)
{
    const size_t arraysize = v.size();
    const It     first     = v.begin();

    if (arraysize == 0) return; // pathological
    size_t i = 0;
    if (howmany >= 4) {
        for (; i <= howmany - N; i += N) {
            It base[N];
            for (int j = 0; j < N; j++)
                base[j] = first;
            size_t n = arraysize;
            while (n > 1) {
                size_t half = n >> 1;
                It     mid[N];
                for (int j = 0; j < N; j++)
                    mid[j] = base[j] + half;
                for (int j = 0; j < N; j++)
                    base[j] = (*mid[j] < targets[i + j]) ? mid[j] : base[j];
                n -= half;
            }
            for (int j = 0; j < N; j++)

                solutions[i + j] = (*base[j] < targets[i + j]) + (base[j] - first);
        }
    }
    for (; i < howmany; i++) {
        binsearch<1>(v, howmany - i, targets + i, solutions + i);
    }
}

template<int N = 4, typename Vec>
static inline void noinline_fun
binsearchpos(Vec& v, int howmany, const size_t* __restrict targets, size_t* __restrict solutions)
{
    const size_t arraysize = v.size();

    if (arraysize == 0) return; // pathological
    size_t i = 0;
    if (howmany >= N) {
        for (; i <= size_t(howmany - N); i += N) {
            size_t base[N];
            for (int j = 0; j < N; j++)
                base[j] = 0;
            size_t n = arraysize;
            while (n > 1) {
                size_t half = n >> 1;
                size_t mid[N];
                for (int j = 0; j < N; j++)
                    mid[j] = base[j] + half;
                for (int j = 0; j < N; j++)
                    base[j] = (v[mid[j]] < targets[i + j]) ? mid[j] : base[j];
                n -= half;
            }
            for (int j = 0; j < N; j++)
                solutions[i + j] = (v[base[j]] < targets[i + j]) + base[j];
        }
    }
    for (; i < howmany; i++) {
        binsearchpos<1>(v, howmany - i, targets + i, solutions + i);
    }
}

template<int N = 4, typename Vec>
static inline void noinline_fun
binsearchpos2(Vec& v, int howmany, const size_t* __restrict targets, size_t* __restrict solutions)
{
    const size_t arraysize = v.size();

    if (arraysize == 0) return; // pathological
    size_t i = 0;
    if (howmany >= N) {
        for (; i <= size_t(howmany - N); i += N) {
            size_t base[N];
            for (int j = 0; j < N; j++)
                base[j] = 0;
            size_t n = arraysize;
            while (n > 1) {
                n >>= 1;

                for (int j = 0; j < N; j++) {
                    size_t b   = base[j];
                    size_t mid = b + n;
                    base[j]    = (v[mid] < targets[i + j]) ? mid : b;
                }
            }
            for (int j = 0; j < N; j++)
                solutions[i + j] = (v[base[j]] < targets[i + j]) + base[j];
        }
    }
    for (; i < howmany; i++) {
        binsearchpos2<1>(v, howmany - i, targets + i, solutions + i);
    }
}

template<int N = 4, typename Vec>
static inline void noinline_fun
binsearchpos3(Vec& v, int howmany, const size_t* __restrict targets, size_t* __restrict solutions)
{
    const size_t        arraysize = v.size();
    const auto          mask      = v.mask;
    const auto          width     = v.width;
    const size_t* const base_ptr  = v.memory.data();

    if (arraysize == 0) return; // pathological
    size_t i = 0;
    if (howmany >= N) {
        for (; i <= size_t(howmany - N); i += N) {
            size_t base[N];
            for (int j = 0; j < N; j++)
                base[j] = 0;

            size_t n = arraysize;
            while (n > 1) {
                n >>= 1;

                for (int j = 0; j < N; j++) {
                    size_t    mid       = base[j] + n;
                    size_t    first_bit = mid * width;
                    bitsize_t idx_bits  = first_bit & 7u;
                    size_t    idx_bytes = first_bit >> 3u;
                    size_t    val_mid;
                    memcpy(&val_mid, reinterpret_cast<const uint8_t*>(base_ptr) + idx_bytes, sizeof(size_t));
                    val_mid = (val_mid >> idx_bits) & mask;
                    base[j] = (val_mid < targets[i + j]) ? mid : base[j];
                }
            }

            for (int j = 0; j < N; j++) {
                solutions[i + j] = (base[j] != 0 || v[0] < targets[i + j]) + base[j];
            }
        }
    }
    for (; i < howmany; i++) {
        binsearchpos3<1>(v, howmany - i, targets + i, solutions + i);
    }
}

template<typename VecArr> using binsearch_t = void (*)(VecArr& v, int, const size_t* __restrict, size_t* __restrict);

template<typename Vec, binsearch_t<Vec> binsearch>
static void
binsearch_bench(benchmark::State& state)
{
    constexpr uint8_t width       = 17;
    constexpr size_t  max_element = (size_t(1) << width) - 1;
    constexpr size_t  log2        = 22;
    constexpr size_t  N           = size_t(1) << log2;

    constexpr size_t query_size = 1024;

    Vec  v(make_sorted_randint_vector(N, width), max_element);
    auto query = make_randint_vector<size_t>(query_size, width);

    for (auto _ : state) {
        std::vector<size_t> solutions(query.size());
        // binsearch4(v.begin(), query.size(), v.size(), query.data(), solutions.data());
        binsearch(v, query.size(), query.data(), solutions.data());
#ifndef NDEBUG
        for (size_t i = 0; i < query.size(); i++) {
            //        std::cout << query[i];
            //        if (solutions[i] >= 1) std::cout << " -" << v[solutions[i] - 1];
            //        if (solutions[i] < v.size())
            //            std::cout << " " << v[solutions[i]];
            //        else
            //            std::cout << " last";
            //        if (solutions[i] + 1 < v.size())
            //            std::cout << " +" << v[solutions[i] + 1];
            //        else
            //            std::cout << " +last";
            //        std::cout << std::endl;

            if (solutions[i] < v.size())
                assert(v[solutions[i]] >= query[i], "wtf");
            else
                assert(v[v.size() - 1] < query[i], "wtf");

            if (solutions[i] >= 1) assert(v[solutions[i] - 1] < query[i], "wtf");
        }
#endif

        benchmark::DoNotOptimize(solutions);
    }
}

BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos<12>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos2<12>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos<12>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos2<12>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos<12>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos2<12>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos3<12>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos<12>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos2<12>);

BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos<8>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos2<8>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos<8>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos2<8>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos<8>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos2<8>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos3<8>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos<8>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos2<8>);

BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos<6>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos2<6>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos<6>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos2<6>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos<6>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos2<6>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos3<6>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos<6>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos2<6>);

BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos<4>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos2<4>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos<4>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos2<4>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos<4>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos2<4>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos3<4>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos<4>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos2<4>);

BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos<1>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_ptrtrick>, binsearchpos2<1>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos<1>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_bitoffset>, binsearchpos2<1>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos<1>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos2<1>);
BENCHMARK_TEMPLATE(binsearch_bench, int_mockup<it_pos>, binsearchpos3<1>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos<1>);
BENCHMARK_TEMPLATE(binsearch_bench, gatbl::int_vector<>, binsearchpos2<1>);

template<typename Vec>
static noinline_fun void
do_iter(Vec& v, size_t mask)
{
    size_t     x   = 0;
    const auto end = v.end();
    for (auto it = v.begin(); it < end; it += 1) {
        x += *it;
        if ((*it & 3) == 0) *it = splitmix64_stateless(*it) & mask;
    }
    benchmark::DoNotOptimize(x);
}

template<typename Vec>
static void
iter(benchmark::State& state)
{
    constexpr size_t log2 = 16;

    constexpr size_t N = size_t(1) << log2;

    uint8_t width = 17;
    use(&width);
    clobber();
    size_t mask = size_t(-1) >> (64 - width);

    Vec v(make_sorted_randint_vector(N, width), mask);

    for (auto _ : state) {
        do_iter(v, mask);
    }
}
BENCHMARK_TEMPLATE(iter, int_mockup<it_ptrtrick>);
BENCHMARK_TEMPLATE(iter, int_mockup<it_bitoffset>);
BENCHMARK_TEMPLATE(iter, int_mockup<it_pos>);
BENCHMARK_TEMPLATE(iter, gatbl::int_vector<>);
