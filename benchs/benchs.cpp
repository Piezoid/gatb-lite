#include <vector>
#include <thread>

#include <immintrin.h>
#include <emmintrin.h>

#include <benchmark/benchmark.h>

#include "gatbl/fastx.hpp"
#include "gatbl/sys/file.hpp"
#include "gatbl/kmer.hpp"

using namespace std;
using namespace gatbl;

using sequence_t = fastq_record<const char*>::substring_t;

template<typename kmer_t = kmer_t, typename Canonical = LexicoCanonical, bool saturated = true>
class kmer_window_valid : private Canonical
{
    using unsaturated = kmer_window_valid<kmer_t, Canonical, false>;

  public:
    using value_type = kmer_t;

    kmer_window_valid(ksize_t k)
      : _mask(~kmer_t(0) >> (CHAR_BIT * sizeof(kmer_t) - 2 * k))
      , _k(k)
      , _left_bitpos(2 * (k - 1))
    {
        clear();
        assume((k & 1) != 0, "k must be odd");
        assume(k <= sizeof(kmer_t) * 4, "k too large");
    }

    kmer_window_valid(unsaturated&& from)
      : _mask(from._mask)
      , _k(from._k)
      , _forward(from._forward)
      , _reverse(from._reverse)
      , _invalid(from._invalid)
    {
        check();
    }

    void clear()
    {
        _reverse = 0;
        _forward = 0;
        _invalid = 0;
    }

    template<typename It>
    hot_fun auto fill(It it) -> decltype(concepts::type_require<It>(concepts::value_require<nuc_t>(*it)))
    {
        clear();

        nuc_t nuc   = *it;
        bool  valid = check_nuc(nuc);

        as_unsaturated().push_back(*it);
        for (ksize_t i = 1; i < _k; ++i) {
            nuc   = *it;
            valid = check_nuc(nuc);
            as_unsaturated().push_back(*it);
        }
        check();
        return it;
    }

    template<typename It>
    hot_fun auto fill(It it)
      -> decltype(concepts::type_require<It>(concepts::value_require<std::pair<nuc_t, bool>>(*it)))
    {
        clear();

        auto nuc_invalid = *it;
        as_unsaturated().push_back(nuc_invalid.first, nuc_invalid.second);
        for (ksize_t i = 1; i < _k; ++i) {
            nuc_invalid = *it;
            as_unsaturated().push_back(nuc_invalid.first, nuc_invalid.second);
        }
        check();
        return it;
    }

    kmer_window_valid& set_forward(kmer_t kmer)
    {
        _forward = kmer;
        _reverse = rcb(kmer, _k);
        _invalid = 0;
        return check();
    }

    kmer_window_valid& set_reverse(kmer_t kmer)
    {
        _reverse = kmer;
        _forward = rcb(kmer, _k);
        _invalid = 0;
        return check();
    }

    static constexpr bool valid_forward = true;
    hot_fun kmer_window_valid& push_back(nuc_t nuc, nuc_t invalid)
    {
        assume(check_nuc(nuc), "invalid nucleotide");
        assume(check_nuc(invalid), "invalid nucleotide");

        _forward <<= 2;
        _forward |= kmer_t(nuc);

        _reverse >>= 2;
        _reverse |= kmer_t(0b10 ^ nuc) << _left_bitpos;

        if (valid_forward) {
            _invalid <<= 2;
            _invalid |= kmer_t(invalid);
        } else {
            _invalid >>= 2;
            _invalid |= kmer_t(invalid) << _left_bitpos;
        }

        if (saturated) {
            _forward &= _mask;
            if (valid_forward) _invalid &= _mask;
        }

        return check();
    }

    hot_fun kmer_window_valid& push_front(nuc_t nuc, nuc_t invalid)
    {
        assume(check_nuc(nuc), "invalid nucleotide");
        assume(check_nuc(invalid), "invalid nucleotide");

        _forward >>= 2;
        _forward |= kmer_t(nuc) << _left_bitpos;

        _reverse <<= 2;
        _reverse |= kmer_t(0b10 ^ nuc);

        if (!valid_forward) {
            _invalid <<= 2;
            _invalid |= kmer_t(invalid);
        } else {
            _invalid >>= 2;
            _invalid |= kmer_t(invalid) << _left_bitpos;
        }

        if (saturated) {
            _reverse &= _mask;
            if (!valid_forward) _invalid &= _mask;
        }

        return check();
    }

    hot_fun kmer_window_valid& push_back(std::pair<nuc_t, bool> nuc_invalid)
    {
        return push_back(nuc_invalid.first, nuc_invalid.second);
    }

    hot_fun kmer_window_valid& push_front(std::pair<nuc_t, bool> nuc_invalid)
    {
        return push_front(nuc_invalid.first, nuc_invalid.second);
    }

    hot_fun kmer_window_valid& push_back(nuc_t nuc)
    {
        bool valid = check_nuc(nuc);
        nuc        = valid ? nuc : 0;
        return push_back(nuc, !valid);
    }

    hot_fun kmer_window_valid& push_front(nuc_t nuc)
    {
        bool valid = check_nuc(nuc);
        nuc        = valid ? nuc : 0;
        return push_front(nuc, !valid);
    }

    ksize_t       size() const { return _k; }
    const kmer_t& forward() const { return _forward; }
    const kmer_t& reverse() const { return _reverse; }
    kmer_t        canon() const { return Canonical::canonize_bidir(_forward, _reverse); }
    bool          valid() const { return _invalid == 0; }

    kmer_t operator*() { return canon(); }

    // private:
    nuc_t check_nuc(nuc_t nuc) const
    {
        assume(nuc < 4 || nuc == nuc_t(-1), "Invalid nuclotide code %u", unsigned(nuc));
        return nuc != nuc_t(-1);
    }

    kmer_window_valid& check()
    {
        assume(_forward <= _mask, "forward kmer is greater than max value");
        assume(_reverse <= _mask, "reverse kmer is greater than max value");
        assume(_invalid <= _mask, "invlid mask is greater than max value");
        assert(!saturated || _forward == rcb(_reverse, _k), "Reversed sequence don't match the forward sequence");
        return *this;
    }

    unsaturated& as_unsaturated() { return reinterpret_cast<unsaturated&>(*this); }

    const kmer_t    _mask;
    kmer_t          _forward, _reverse, _invalid;
    const ksize_t   _k;
    const bitsize_t _left_bitpos;
};

std::vector<sequence_t>
read_sample(const char* filename = "/home/mkerbiri/data/frag_1.fastq", size_t n = 4096)
{
    sys::file_descriptor fd(filename);
    auto                 data = fd.mmap<const char>();
    // FIXME: those are optional but we really want to knwow if they work for testing
    data.advise_hugepage();
    data.advise_sequential();

    using range_t = sequence_range<fastq_record<const char*>>;
    range_t range = data;

    std::vector<sequence_t> sequences;
    sequences.reserve(n);
    RANGES_FOR(auto& rec, range)
    {
        if (n-- == 0) break;
        sequences.emplace_back(rec.sequence());
    }
    return sequences;
}

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

static void
Baseline(benchmark::State& state)
{
    auto samples = read_sample();
    for (auto _ : state) {
        for (auto& x : samples)
            use(x);
    }
}
BENCHMARK(Baseline);

static void
KmerIter1(benchmark::State& state)
{
    auto samples = read_sample();

    ksize_t k = 31;
    use(k);
    clobber();
    kmer_window_valid<> win(k);
    for (auto _ : state) {
        size_t hash = 0;
        for (auto& rec : samples) {
            using seq_t = remove_reference_t<decltype(rec)>;
            dna_ascii_range<seq_t> seq(rec);
            for (auto kmer : seq | win)
                hash += win.valid() + win.canon();
        }
        use(hash);
    }
}
BENCHMARK(KmerIter1);

__attribute__((always_inline)) inline __m128i
convert(__m128i& __restrict in_out)
{
    // Take an upper-case snapshot of input
    auto inU = _mm_and_si128(in_out, _mm_set1_epi8(uint8_t(223)));
    // Convert to 2bit: (x >> 1) & 0b11
    in_out = _mm_and_si128(_mm_srli_epi64(in_out, 1), _mm_set1_epi8(uint8_t(3)));

    // Send 'A's to 'C', then test C
    auto isAorC = _mm_cmpeq_epi8(_mm_or_si128(inU, _mm_set1_epi8(2)), _mm_set1_epi8('C'));
    auto isT    = _mm_cmpeq_epi8(inU, _mm_set1_epi8('T')); // Test 'T'
    auto isG    = _mm_cmpeq_epi8(inU, _mm_set1_epi8('G')); // Test 'G'
    // Or the three results
    return _mm_or_si128(_mm_or_si128(isT, isG), isAorC);
}

__attribute__((always_inline)) inline uint16_t
convert(__m128i& out, const char* input, size_t len)
{
    __m128i in;
    if (len == 16) {
        memcpy(&in, input, 16);
    } else {
        memset(&in, 0, sizeof(__m128i));
        for (unsigned i = 0; i < len; i++) {
            reinterpret_cast<char*>(&in)[i] = input[i];
        }
    }

    auto valid = convert(in);
    out        = in;
    return uint16_t(~_mm_movemask_epi8(valid));
}

static void
NucIter1(benchmark::State& state)
{
    auto samples = read_sample();
    for (auto _ : state) {
        for (auto& x : samples) {

            auto*   it  = x.begin();
            auto*   end = x.end();
            __m128i out;
            while (true) {
                if (it + 16 <= end) {
                    size_t invalid = convert(out, it, 16);
                    use(out);
                    use(invalid);
                    it += 16;
                } else {
                    size_t invalid = convert(out, it, end - it);
                    use(out);
                    use(invalid);
                    break;
                }
            }
        }
    }
}
BENCHMARK(NucIter1);

struct nuc_iter
{
    using char_iterator_traits = std::iterator_traits<const char*>;
    using value_type           = std::pair<nuc_t, bool>;
    using reference            = value_type&;
    using iterator_category    = typename char_iterator_traits::iterator_category;
    using difference_type      = typename char_iterator_traits::difference_type;
    using pointer              = typename char_iterator_traits::pointer;

    using vec_t                        = __m128i;
    static constexpr size_t vec_size   = sizeof(vec_t);
    static constexpr size_t vec_n      = 2;
    static constexpr size_t batch_size = vec_n * vec_size;

    nuc_iter(const char* it, const char* end)
      : _it(it)
      , _end(end)
      , _pos(0)
    {
        assume(_it + batch_size <= _end, "Too small input");
        next();
    }
    nuc_iter() = default;

    forceinline_fun void next()
    {
        size_t size = _end - _it;
        assume(size > 0, "Zero size next()");
        _pos = size >= batch_size ? 0 : batch_size - size;
        _it -= _pos;

        _invalid = 0;
        for (unsigned i = 0; i < vec_n; i++) {
            _invalid <<= vec_size;
            _invalid |= convert(_nucs[i], _it, vec_size);
            _it += vec_size;
        }
        _invalid >>= _pos;
    }

    forceinline_fun nuc_iter& operator++()
    {
        _invalid >>= 1;
        if (++_pos >= batch_size && _it != _end) next();
        return *this;
    }

    std::pair<nuc_t, bool> operator*()
    {
        assume(_pos < batch_size, "_pos out of range : %u", _pos);
        return {reinterpret_cast<uint8_t(&)[batch_size]>(_nucs[0])[_pos], _invalid & 1};
    }

    bool operator!=(const char* end)
    {
        assume(end == _end, "Testing different end");
        return _pos < batch_size;
    }

    __m128i       _nucs[vec_n];
    const char*   _it;
    const char*   _end;
    uint_fast64_t _invalid;
    uint_fast16_t _pos;
};

static void
NucIter2(benchmark::State& state)
{
    auto samples = read_sample();

    nuc_iter iter;
    for (auto _ : state) {
        for (auto& seq : samples) {
            auto* end = seq.end();
            // std::cout << std::string(seq.begin(), seq.end()) << endl;
            for (iter = nuc_iter(seq.begin(), end); iter != end; ++iter) {
                auto  p       = *iter;
                nuc_t nuc     = p.first;
                bool  invalid = p.second;
                // std::cerr << (invalid ? 'N' : "ACTG"[nuc]);
                use(nuc);
                use(invalid);
            }
            // std::cout << endl;
        }
    }
}
BENCHMARK(NucIter2);

static void
KmerIter2(benchmark::State& state)
{
    auto    samples = read_sample();
    ksize_t k       = 7;
    use(k);
    clobber();
    kmer_window_valid<uint32_t> win(k);

    for (auto _ : state) {
        size_t hash = 0;
        for (auto& rec : samples) {
            iterator_pair<nuc_iter, const char*> seq{{rec.begin(), rec.end()}, rec.end()};
            for (auto kmer : seq | win)
                hash += win.valid() + win.canon();
        }
        use(hash);
    }
}
BENCHMARK(KmerIter2);

size_t
aligned_right(void*& __ptr, size_t align)
{
    const auto __intptr  = reinterpret_cast<uintptr_t>(__ptr);
    const auto __aligned = (__intptr - 1u + align) & -align;
    __ptr                = reinterpret_cast<void*>(__aligned);
    return __aligned - __intptr;
}

size_t
aligned_left(void*& __ptr, size_t align)
{
    const auto __intptr  = reinterpret_cast<uintptr_t>(__ptr);
    const auto __aligned = __intptr & -align;
    __ptr                = reinterpret_cast<void*>(__aligned);
    return __intptr - __aligned;
}

static void flatten_fun
            KmerIter3(benchmark::State& state)
{
    auto  samples = read_sample();
    auto* scratch = reinterpret_cast<__m128i*>(aligned_alloc(16, 4096));

    ksize_t k = 7;
    use(k);
    clobber();

    assume(k <= 16, "wut");

    for (auto _ : state) {
        size_t hash = 0;
        for (auto& rec : samples) {
            kmer_window_valid<uint32_t, LexicoCanonical> win(k);
            auto                                         it  = rec.begin();
            const auto                                   end = rec.end();

            __m128i _v      = _mm_lddqu_si128(reinterpret_cast<const __m128i*>(it));
            auto    invalid = (~reinterpret_cast<__v16qu>(convert(_v))) & 3;
            auto    v       = reinterpret_cast<__v16qu>(_v);

            //            if (memchr(it, 'N', end - it) != nullptr) { // He
            //                cout << "blih" << endl;
            //            }

            for (int j = 0; j < k; j++) {
                win.as_unsaturated().push_back(v[j], invalid[j]);
                //                if (invalid[k]) cerr << "blah";
            }
            win.check();

            hash += win.valid() + win.canon();

            it += k;
            size_t len = end - it;
            memcpy(scratch, it, len);
            memset(reinterpret_cast<char*>(scratch) + len, 0, 16);

            len = (len + (16 - 1)) / 16;
            for (unsigned i = 0; i < len; i++) {
                __m128i _v = scratch[i];
                invalid    = (~reinterpret_cast<__v16qu>(convert(_v))) & 3;
                v          = reinterpret_cast<__v16qu>(_v);

                //#pragma unroll 16
                for (int j = 0; j < 16; j++) {

                    //                    if (invalid[j]) { // ho
                    //                        cerr << "blah";
                    //                    }
                    win.push_back(v[j], invalid[j]);
                    hash += win.valid() + win.canon();

                    //                    v       = _mm_bsrli_si128(v, 1);
                    //                    invalid = _mm_bsrli_si128(invalid, 1);
                }
            }
        }
        use(hash);
    }

    free(scratch);
}
BENCHMARK(KmerIter3);

static void flatten_fun
            KmerIter4(benchmark::State& state)
{
    auto  samples = read_sample();
    auto* scratch = reinterpret_cast<__m128i*>(aligned_alloc(16, 4096));

    ksize_t k = 7;
    use(k);
    clobber();

    assume(k <= 16, "wut");

    for (auto _ : state) {
        size_t hash = 0;
        for (auto& rec : samples) {
            kmer_window_valid<uint32_t, LexicoCanonical> win(k);
            auto                                         it  = rec.begin();
            const auto                                   end = rec.end();

            __m128i v       = _mm_lddqu_si128(reinterpret_cast<const __m128i*>(it));
            __m128i invalid = _mm_andnot_si128(convert(v), _mm_set1_epi8(3));

            //            if (memchr(it, 'N', end - it) != nullptr) { // He
            //                cout << "blih" << endl;
            //            }

            for (int j = 0; j < k; j++) {
                win.as_unsaturated().push_back(_mm_extract_epi8(v, 0), _mm_extract_epi8(invalid, 0));
                v       = _mm_bsrli_si128(v, 1);
                invalid = _mm_bsrli_si128(invalid, 1);
                //                if (invalid[k]) cerr << "blah";
            }
            win.check();

            hash += win.valid() + win.canon();

            it += k;
            size_t len = end - it;
            memcpy(scratch, it, len);
            memset(reinterpret_cast<char*>(scratch) + len, 0, 16);

            len = (len + (16 - 1)) / 16;
            for (unsigned i = 0; i < len; i++) {
                v       = scratch[i];
                invalid = _mm_andnot_si128(convert(v), _mm_set1_epi8(3));

#pragma unroll 16
                for (int j = 0; j < 16; j++) {

                    //                    if (invalid[j]) { // ho
                    //                        cerr << "blah";
                    //                    }
                    win.push_back(_mm_extract_epi8(v, 0), _mm_extract_epi8(invalid, 0));
                    hash += win.valid() + win.canon();

                    v       = _mm_bsrli_si128(v, 1);
                    invalid = _mm_bsrli_si128(invalid, 1);
                }
            }
        }
        use(hash);
    }

    free(scratch);
}
BENCHMARK(KmerIter4);
