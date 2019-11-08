#include <vector>
#include <thread>

#include <sys/mman.h>

#include <immintrin.h>
#include <emmintrin.h>

#include <benchmark/benchmark.h>

#include "gatbl/fastx.hpp"
#include "gatbl/sys/file.hpp"
#include "gatbl/kmer.hpp"

#undef assert
#include <cassert>

using namespace std;
using namespace gatbl;

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

// Westmere Xeon X8675
static constexpr size_t line_size = 64;
static constexpr size_t l1_size   = size_t(32) << 10;
static constexpr size_t l2_size   = size_t(256) << 10;
static constexpr size_t l3_size   = size_t(8) << 20;

struct CacheLevel
{
    using vec_t                         = __m128i;
    static constexpr size_t align       = l1_size; // >64
    static constexpr size_t nelems_line = line_size / sizeof(vec_t);

    static_assert(line_size % sizeof(vec_t) == 0, "fractional number of words in a cache line");
    static_assert(align >= line_size, "sub-aligned to cache lines");

    /// Allocate a blob of `size` bytes
    CacheLevel(size_t size)
      : _size(size)
    {
        void* vptr;
        _size &= size_t(-align);
        posix_memalign(&vptr, align, _size);
        if (_size >= 2 << 20) assert(posix_madvise(vptr, _size, MADV_HUGEPAGE) == 0);
        _ptr = reinterpret_cast<vec_t*>(vptr);
        check_line_aligned(_ptr);
    }

    void memset(uint8_t val = 0xaa, size_t nlines_batch = 1)
    {
        const size_t nelems_batch = nlines_batch * nelems_line;

        const vec_t vec = _mm_set1_epi8(val);

        auto* it   = begin();
        auto* last = end();
        for (; it < last;) {
            check_line_aligned(it);
            auto* last_inbatch = it + nelems_batch;
            for (; it < last_inbatch; ++it) {
                _mm_store_si128(it, vec);
            }
        }

        clobber();

        clobber();
    }

    void memsetnt(uint8_t val = 0xaa, size_t nlines_batch = 1)
    {
        const size_t nelems_batch = nlines_batch * nelems_line;

        const vec_t vec = _mm_set1_epi8(val);

        auto* it   = begin();
        auto* last = end();
        for (; it < last;) {
            check_line_aligned(it);
            auto* last_inbatch = it + nelems_batch;
            for (; it < last_inbatch; ++it) {
                _mm_stream_si128(it, vec);
            }
        }

        clobber();
    }

    void clflush()
    {
        auto* last = end();
        for (auto* it = begin(); it < last; it += line_size) {
            check_line_aligned(it);
            _mm_clflush(it);
        }
        clobber();
    }

    void load(size_t nlines_batch = 1)
    {
        const size_t nelems_batch = nlines_batch * nelems_line;

        auto* it   = begin();
        auto* last = end();
        vec_t hash = {0, 0};

        for (; it < last;) {
            check_line_aligned(it);
            auto* last_inbatch = it + nelems_batch;
            for (; it < last_inbatch; ++it) {
                hash += _mm_load_si128(it);
            }
        }

        use_vec(hash);
    }

    void loadnt(size_t nlines_batch = 1)
    {
        const size_t nelems_batch = nlines_batch * nelems_line;

        auto* it   = begin();
        auto* last = end();
        vec_t hash = {0, 0};

        for (; it < last;) {
            check_line_aligned(it);
            auto* last_inbatch = it + nelems_batch;
            for (; it < last_inbatch; ++it) {
                hash += _mm_stream_load_si128(it);
            }
        }

        use_vec(hash);
    }

    void load_prefetchnt(size_t nlines_batch = 1, size_t nlines_lookahead = 32)
    {
        const size_t nelems_batch     = nlines_batch * nelems_line;
        const size_t nelems_lookahead = nlines_lookahead * nelems_batch;

        auto* it   = begin();
        auto* last = end();
        vec_t hash = {0, 0};

        for (; it < last;) {
            auto* last_inbatch = it + nelems_batch;
            for (auto* it_inbatch = it; it_inbatch < last_inbatch; it_inbatch += nelems_line) {
                check_line_aligned(it);
                _mm_prefetch(nelems_lookahead + it_inbatch, _MM_HINT_NTA);
            }

            for (; it < last_inbatch; ++it) {
                hash += _mm_load_si128(it);
            }
        }

        use_vec(hash);
    }

    ~CacheLevel() { free(_ptr); }

    size_t size() const { return _size & size_t(-align); }
    vec_t* begin() const { return reinterpret_cast<vec_t*>(reinterpret_cast<uintptr_t>(_ptr) & size_t(-align)); }
    vec_t* end() const { return begin() + size() / sizeof(vec_t); }

  private:
    void use_vec(vec_t& hash)
    {
        int64_t hash0 = _mm_extract_epi64(hash, 0);
        int64_t hash1 = _mm_extract_epi64(hash, 1);
        ::use(hash0);
        ::use(hash1);
    }

    void clobber()
    {
        ::use(_ptr);
        ::clobber();
    }

    template<typename T> static void check_line_aligned(const T* ptr)
    {
        assert(reinterpret_cast<uintptr_t>(ptr) % line_size == 0);
    }

    vec_t* _ptr;
    size_t _size;
};

constexpr size_t stream_size = l3_size;
constexpr size_t cache_size  = l2_size;

constexpr size_t nlines_batch     = 2;
constexpr size_t nlines_lookahead = 8;

static void
load_cache(benchmark::State& state)
{
    CacheLevel stream(stream_size);
    stream.memsetnt(0x55);
    stream.clflush();

    CacheLevel cache(cache_size);
    cache.memset(0xaa);

    for (auto _ : state) {
        state.PauseTiming();
        clobber();
        state.ResumeTiming();

        cache.load(nlines_batch);
    }
}
BENCHMARK(load_cache);

static void
flush_then_load_cache(benchmark::State& state)
{
    CacheLevel stream(stream_size);
    stream.memsetnt(0x55);
    stream.clflush();

    CacheLevel cache(cache_size);
    cache.memset(0xaa);

    for (auto _ : state) {
        state.PauseTiming();
        cache.clflush();
        state.ResumeTiming();

        cache.load(nlines_batch);
    }
}
BENCHMARK(flush_then_load_cache);

static void
load_stream_then_load_cache(benchmark::State& state)
{
    CacheLevel stream(stream_size);
    stream.memsetnt(0x55);
    stream.clflush();

    CacheLevel cache(cache_size);
    cache.memset(0xaa);

    for (auto _ : state) {
        state.PauseTiming();
        stream.load(nlines_batch);
        state.ResumeTiming();

        cache.load(nlines_batch);
    }
}
BENCHMARK(load_stream_then_load_cache);

static void
loadnt_stream_then_load_cache(benchmark::State& state)
{
    CacheLevel stream(stream_size);
    stream.memsetnt(0x55);
    stream.clflush();

    CacheLevel cache(cache_size);
    cache.memset(0xaa);

    for (auto _ : state) {
        state.PauseTiming();
        stream.loadnt(nlines_batch);
        state.ResumeTiming();

        cache.load(nlines_batch);
    }
}
BENCHMARK(loadnt_stream_then_load_cache);

static void
loadprefetchnt_stream_then_load_cache(benchmark::State& state)
{
    CacheLevel stream(stream_size);
    stream.memsetnt(0x55);
    stream.clflush();

    CacheLevel cache(cache_size);
    cache.memset(0xaa);

    for (auto _ : state) {
        state.PauseTiming();
        stream.load_prefetchnt(nlines_batch, nlines_lookahead);
        state.ResumeTiming();

        cache.load(nlines_batch);
    }
}
BENCHMARK(loadprefetchnt_stream_then_load_cache);

static void
memset_stream_then_load_cache(benchmark::State& state)
{
    CacheLevel stream(stream_size);
    stream.memsetnt(0x55);
    stream.clflush();

    CacheLevel cache(cache_size);
    cache.memset(0xaa);

    for (auto _ : state) {
        state.PauseTiming();
        stream.memset(0x55, nlines_batch);
        state.ResumeTiming();

        cache.load(nlines_batch);
    }
}
BENCHMARK(memset_stream_then_load_cache);

static void
memsetnt_stream_then_load_cache(benchmark::State& state)
{
    CacheLevel stream(stream_size);
    stream.memsetnt(0x55);
    stream.clflush();

    CacheLevel cache(cache_size);
    cache.memset(0xaa);

    for (auto _ : state) {
        state.PauseTiming();
        stream.memsetnt(0x5a, nlines_batch);
        state.ResumeTiming();

        cache.load(nlines_batch);
    }
}
BENCHMARK(memsetnt_stream_then_load_cache);
