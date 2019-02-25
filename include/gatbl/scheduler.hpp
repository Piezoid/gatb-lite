#ifndef SCHEDULER_HPP
#define SCHEDULER_HPP

#include <cstring>
#include <memory>
#include <algorithm>
#include <functional>

#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "gatbl/utils/compatibility.hpp"
#include "gatbl/common.hpp"
#include "gatbl/sys/thread.hpp"

namespace gatbl {
namespace memory {

/// Wrapper of a T that might be allocated on heap if it is larger than max_size
template<typename T, size_t max_size, bool remote = (sizeof(T) > max_size)> struct bounded_size_wrapper;

template<typename T, size_t max_size> struct bounded_size_wrapper<T, max_size, false>
{
    using value_type                       = T;
    static constexpr bool have_indirection = false;

    template<typename... Args>
    bounded_size_wrapper(Args&&... args) noexcept(noexcept(value_type(std::forward<Args>(args)...)))
      : value(std::forward<Args>(args)...)
    {}
    bounded_size_wrapper(bounded_size_wrapper&&) = default;
    value_type&       get() noexcept { return value; }
    const value_type& get() const noexcept { return value; }

  private:
    value_type value;
};

template<typename T, size_t max_size> struct bounded_size_wrapper<T, max_size, true>
{
    using value_type                       = T;
    static constexpr bool have_indirection = true;

    template<typename... Args>
    bounded_size_wrapper(Args&&... args) noexcept(noexcept(value_type(std::forward<Args>(args)...)))
      : value(new (std::nothrow) T(std::forward<Args>(args)...))
    {}
    value_type&       get() noexcept { return *value; }
    const value_type& get() const noexcept { return *value; }
    ~bounded_size_wrapper() noexcept(noexcept(std::declval<value_type*>->~value_type())) { delete value; }

  private:
    value_type* value;
};

} // namespace memory

namespace sched {

constexpr size_t hardware_destructive_interference_size = 64;

template<typename Fun, size_t cache_lines = 1> class jobslot;

template<typename... Args, size_t cache_lines>
class alignas(hardware_destructive_interference_size) jobslot<void(Args...), cache_lines>
{
  public:
    class jobtoken;

  protected:
    using func_ptr           = void (*)(jobtoken&&, Args...) CPP17_NOEXCEPT;
    using finish_handler_ptr = void (*)(void*) CPP17_NOEXCEPT;

  private:
    static constexpr size_t data_size
      = cache_lines * hardware_destructive_interference_size - sizeof(std::atomic<size_t>) - sizeof(func_ptr);
    template<typename F> using data_wrapper = memory::bounded_size_wrapper<F, data_size>;

    alignas(hardware_destructive_interference_size) char _data[data_size];

    // The triaval constructor make an empty job in finished state
    std::atomic<size_t> _state_idx{0};
    // _state_idx = 0 is empty job ready to be allocated in the pool
    // _state_idx = 1 is a job being deleted
    // _state_idx = 1 + n is scheduled job with (n-1) alive references
    func_ptr _function = nullptr;

  public:
    jobslot() noexcept      = default;
    jobslot(const jobslot&) = delete;
    jobslot& operator=(const jobslot&) = delete;

#if DEBUG
    ~jobslot() noexcept { assert(empty(), "non executed job destroyed"); }
#else
    ~jobslot() noexcept = default;
#endif

    bool empty(std::memory_order mo = std::memory_order_acquire) const noexcept
    {
        bool isempty = _state_idx.load(mo) == 0;
        assert(!isempty || _function == nullptr, "Disowned and ran job with a live function pointer");
        return isempty;
    }

    class jobtoken
    {
      private:
        friend class jobslot;
        jobslot* _j = nullptr;

        void set_noacquire(jobslot* j) noexcept { _j = j; }

        const jobslot& slot() const noexcept
        {
            assume(operator bool(), "Accessing a void job token");
            return *_j;
        }

        jobslot& slot() noexcept
        {
            assume(operator bool(), "Accessing a void job token");
            return *_j;
        }

      protected:
        void run(Args... args) noexcept
        {
            jobslot& s = slot();
            // Precondition : have a function defined and alive jobrefs
            assert(s._state_idx.load(std::memory_order_acquire) >= 2, "Trying to run a disowned job");
            assert(*this, "Trying to run an empty job");
            s._function(std::move(*this), std::forward<Args>(args)...);
        }

        void init(func_ptr f) { slot()._function = f; }

      public:
        /// An empty jobref not referencing any jobconst
        jobtoken() noexcept = default;

        /// Duplicate ownership of the job
        explicit jobtoken(const jobtoken& other) noexcept
          : _j(other._j)
        {
            if (operator bool()) {
                size_t state_idx = _j->_state_idx.fetch_add(1, std::memory_order_relaxed);
                assume(state_idx > 0, "acquire() on unititialized job");
            }
        }

        /// Pass on ownership of the job
        jobtoken(jobtoken&& other) noexcept
          : _j(other._j)
        {
            other._j = nullptr;
        }

        jobtoken& operator=(jobtoken&& other)
        {
            std::swap(_j, other._j);
            return *this;
        }

        /// Release the ownership of the job
        ~jobtoken() noexcept
        {
            if (operator bool()) {
                // Decrement the sate counter
                auto state_idx = _j->_state_idx.fetch_sub(1,
                                                          // TSan doesn't handle well the conditional fence here
                                                          __has_feature(thread_sanitizer) ? std::memory_order_acq_rel
                                                                                          : std::memory_order_release);

                assume(state_idx > 1, "Double free job slot");
                if (state_idx == 2) // If it reaches 1 the job completion handler is executed
                {
                    std::atomic_thread_fence(std::memory_order_acquire);
                    // Publishing _state_idx == 1 allows to postpone the reuse of the pool's slot while running the
                    // completion chain.
                    if (_j->_function != nullptr) {
                        // Call the completion handler
                        reinterpret_cast<finish_handler_ptr>(_j->_function)(_j->_data);
                        if (DEBUG) _j->_function = nullptr;
                    }

                    // Mark the job as finished so it can be reallocated.
                    // Pretty sure that relaxed ordering is enough here:
                    // We don't mind if the allocator sees the previously released value _state_idx=1 instead of 0.
                    // Sadly, TSan thinks otherwise...
                    static constexpr std::memory_order mo
                      = __has_feature(thread_sanitizer) ? std::memory_order_release : std::memory_order_relaxed;
                    if (DEBUG) {
                        // printf("%lu finished job %p\n", std::this_thread::get_id(), this);
                        assert(_j->_state_idx.compare_exchange_strong(--state_idx, 0, mo, std::memory_order_relaxed),
                               "Race condition in job release()");
                    } else {
                        _j->_state_idx.store(0, mo);
                    }
                }
                _j = nullptr;
            }
        }

        /// Returns true the instance is referencing a job
        operator bool() const noexcept
        {
            if (_j != nullptr) {
                assert(not _j->empty(), "Empty job referenced");
                return true;
            } else
                return false;
        }

        size_t owners_count(std::memory_order mo = std::memory_order_acquire) const noexcept
        {
            size_t state_idx = slot()._state_idx.load(mo);
            // There should be at least one owner (including this instance):
            assert(state_idx >= 2, "Acessing a disowned job through a owning ref");
            return state_idx - 1;
        }

        bool single_owner(std::memory_order mo = std::memory_order_acquire) const noexcept
        {
            return owners_count(mo) == 1;
        }
    };

    // Typed interface over a job token: allows to access and mutate the slot content as a Functor instance
    template<typename F, typename Base = jobtoken> struct typed : public Base
    {
        static constexpr inline bool have_indirection = data_wrapper<F>::have_indirection;

        F& functor() { return reinterpret_cast<data_wrapper<F>*>(Base::slot()._data)->get(); }

        const F& functor() const { return reinterpret_cast<data_wrapper<F>*>(Base::slot()._data)->get(); }

        using Base::Base;

      protected:
        friend class jobslot;

        /// Downcasting from an untyped token: as dangerous as reinterpret_cast !
        template<typename... _Args>
        typed(_Args&&... args)
          : Base(std::forward<_Args>(args)...)
        {}

        template<typename... ConstrArgs> void emplace(ConstrArgs&&... constr_args)
        {
            new (Base::slot()._data) data_wrapper<F>(std::forward<ConstrArgs>(constr_args)...);
        }

        /// Set deletion handler (should only be used once the job is startedworker_pool (e.g. at epilogue of the job
        /// function)
        void schedule_destructor()
        {
            finish_handler_ptr finish = std::is_trivially_destructible<F>::value ? nullptr : [](void* ptr) noexcept
            {
                reinterpret_cast<F*>(ptr)->~F();
            };
            Base::slot()._function = reinterpret_cast<func_ptr>(finish);
        }

        /// An example of job initialization that can be achieved with the above mutator
        template<typename... ConstrArgs> void init(ConstrArgs&&... args)
        {
            Base::init([](jobtoken && token, Args... args) noexcept {
                // In the function, the type of the functor is forgoten, so get it back by constructing a typed token of
                // the right type
                typed tytk(std::move(token));
                auto& functor = tytk.functor();
                tytk.schedule_destructor();
                functor(std::move(tytk));
            });
            this->emplace(std::forward<ConstrArgs>(args)...);
        }
    };

    /// Initialize the jobslot and return multiple job_ptr in one go
    template<size_t n, typename Token, typename... ConstrArgs> std::array<Token, n> init(ConstrArgs&&... constr_args)
    {

        static_assert(n > 0, "At least one job reference must be instanciated");

        // jobslot type invariants are checked here, now that the type is complete
        static_assert(sizeof(jobslot) == cache_lines * hardware_destructive_interference_size,
                      "job is the size of some cache line");
        static_assert(alignof(jobslot) % hardware_destructive_interference_size == 0,
                      "job is aligned on cache line boundaries");
        assert(reinterpret_cast<uintptr_t>(this) % hardware_destructive_interference_size == 0,
               "Unaligned job slot"); // RT version of the above

        // State transition to an occupied jobslot
        assert(empty(), "Allocated a non freed job");
        _state_idx.store(1 + n, std::memory_order_relaxed);

        std::array<Token, n> refs{};
        for (Token& ref : refs)
            ref.set_noacquire(this);

        refs[0].init(std::forward<ConstrArgs>(constr_args)...);

        return refs;
    }
};

/** Concurrent SPMC stack/queue for work stealing
 *
 * The producer-owner can push and pop in FIFO order from one end while other thread might steal elements in LIFO order
 * on the other end. _top and _bottom atomic indexes are read and write positions of a ring buffer exploiting the
 * modular arithmetic of the unsigned indice type.
 * @see https://www.snellman.net/blog/archive/2016-12-13-ring-buffers/
 *
 *
 */
template<typename T, size_t cap_bits> struct work_stealing_queue
{
    using idx_t                     = size_t;
    static constexpr idx_t capacity = 1ULL << cap_bits;
    static constexpr idx_t mask     = capacity - 1;

    /** Try to push a item (reserved to the single-producer/owner)
     * In case of success (queue is not full) the ownership of the item is passed to the queue and
     * a trivially constructed element is returned. Otherwise the element is returned to the caller.
     */
    T push(T&& job) noexcept
    {
        assert(job, "pushing empty job");
        // Bottom: relaxed load, since only _bottom is only modified by the worker owning the queue
        idx_t b = _bottom.load(std::memory_order_relaxed);

        // Top: relaxed load, since _top is monotonically increasing we only risk false
        // positives in the "is full" test bellow
        idx_t t = _top.load(std::memory_order_relaxed);

        if (likely(b - t < capacity)) { // not fullatomic_thread_fence
            _arr[b++ & mask] = std::move(job);
            _bottom.store(b, std::memory_order_release);
            // printf("%lu pushed job %p in slot %lu, top=%lu bottom=%lu=>%lu\n", std::this_thread::get_id(), &job, (b -
            // 1) & mask, t, b - 1, b);
            return {};
        }
        return std::move(job);
    }

    /// Pops in FIFO order (reserved to the single-producer/owner)
    T pop() noexcept
    {
        // First we decrement _bottom: load is relaxed since we are the only thread to update it,
        // but publishes the updated value for stealers
        idx_t bottom = _bottom.fetch_sub(1, std::memory_order_release);

        idx_t top = _top.load(std::memory_order_acquire);
        //        printf("%lu pop state b:%lu t:%lu\n", id, bottom, t0);
        T item = {};

        if (likely(bottom - top > 1)) { // More than one job left
            item = std::move(_arr[(bottom - 1) & mask]);
            // printf("%lu get job %p at slot %lu (top=%lu)\n", id, item, (bottom - 1) & mask, top);
            assume(item, "poping invalid item");
            return item;
        }

        if (bottom - top
            == 1) { // Last job in the queue: it could also be extracted concurrently from the other end by steal()
            // So we also do a kind of steal() here, increasing top with CMPXCHG
            if (_top.compare_exchange_weak(top, top + 1, std::memory_order_acq_rel, std::memory_order_relaxed)) {
                // we were alone incrementing top
                item = std::move(_arr[(bottom - 1) & mask]);
                assume(item, "poping invalid item");
                // printf("%lu get job %p at slot %lu (top=%lu)\n", id, item, (bottom - 1) & mask, top);
            }
        }

        // restore bottom
        _bottom.store(bottom, std::memory_order_release);
        assert(_bottom == _top, "Unexpected top idx change");
        return item;
    }

    /// Steal an element in LIFO order. Can be called from any thread.
    T steal() noexcept
    {
        // Relaxed load since the value will be checked with a CAS
        idx_t t = _top.load(std::memory_order_relaxed);
        idx_t b = _bottom.load(std::memory_order_acquire);

        //        printf("%lu steal state b:%lu t:%lu\n", id, b, t0);

        // The queue not empty if bottom > top. Howerver, we can't test that because we want to tolerate wrapparound of
        // the counter empty => bottom == top, unless bottom is temporarly decremented by the queue owner inspecting it
        // with pop() in wich case t == b + 1.
        if (t != b && t != b + 1) {
            // the compare_exchange_weak function serves as a compiler barrier, and guarantees that the read happens
            // before the CAS.
            if (_top.compare_exchange_weak(t, t + 1, std::memory_order_acq_rel, std::memory_order_relaxed)) {
                return std::move(_arr[t & mask]);
            }
        }

        // empty queue
        return {};
    }

    int approximate_size() const noexcept
    { // TODO: check if its a win to avoid a memory fence here (might loose oportunities for stealing)
        return int(_bottom.load(std::memory_order_relaxed)) - int(_top.load(std::memory_order_relaxed));
    }

#if DEBUG
    ~work_stealing_queue() { assert(_bottom == _top, "Jobs remaining in the queue upon destruction"); }
#endif

  private:
    alignas(hardware_destructive_interference_size) T _arr[capacity];
    alignas(hardware_destructive_interference_size) std::atomic<idx_t> _bottom{0};
    alignas(hardware_destructive_interference_size) std::atomic<idx_t> _top{0};
};

/**
 * A simple allocation pool with inline storage for objects with an empty state.
 * Allocation is done by scanning the pool starting after the last allocated slot.
 *
 * The (non-standard) allocate() function (without args) may give up and return a nullptr.
 */
template<typename T, size_t bits = 6> struct scanning_pool_allocator
{
    using value_type             = T;
    using pointer                = value_type*;
    static constexpr size_t size = 1ULL << bits;

    constexpr size_t max_size() const noexcept { return size; }

    /// Tries to allocate one object in the first found slot
    value_type* allocate() noexcept
    {
        assert((uintptr_t)(_pool) % hardware_destructive_interference_size == 0, "Unaligned");
        pointer slot0 = &_pool[_idx++ & _mask];
        pointer slot  = slot0;
        bool    empty = slot->empty();
        while (unlikely(!empty)) {
            slot  = &_pool[_idx++ & _mask];
            empty = slot->empty();
            if (unlikely(!empty && slot == slot0)) // We looped back to the first tried slot, give up for now.
                return nullptr;
        }
        return slot;
    }

    void allocate(size_t n)
    {
        if (n == 1) throw std::bad_alloc();
        pointer p = allocate();
        if (p == nullptr) throw std::bad_alloc();
    }

    bool all_empty() const noexcept
    {
        bool res = true;
        for (const value_type& v : _pool)
            res &= v.empty();
        return res;
    }

    CPP14_CONSTEXPR scanning_pool_allocator() noexcept { assert(all_empty(), "Non empty pool at construction"); }

    ~scanning_pool_allocator() { assert(all_empty(), "Non empty pool at destruction"); }

  private:
    static constexpr size_t _mask = size - 1;
    value_type              _pool[size];
    size_t                  _idx = 0;
};

template<size_t cap_bits = 6, size_t pool_bits = cap_bits> class wker
{
    using jobslot                                        = jobslot<void(wker&)>;
    using base_jobtoken                                  = typename jobslot::jobtoken;
    template<typename F, typename Base> using base_typed = typename jobslot::template typed<F, Base>;

  public:
    struct ctx;
    struct jobtoken : public base_jobtoken
    {
        using Base = typename jobslot::jobtoken;
        using base_jobtoken::base_jobtoken;

      protected:
        template<typename Tk>
        jobtoken(Tk&& from)
          : base_jobtoken(std::forward<Tk>(from)){};

        friend class wker;
        friend jobslot;
    };

    template<typename F, typename Base = jobtoken> struct typed : public base_typed<F, Base>
    {
        using base_t = base_typed<F, Base>;
        using base_t::base_t;

        // Return a new reference to the job
        typed<F, jobtoken> token() const { return typed<F, jobtoken>(*this); }

      protected:
        friend class wker;
        friend jobslot;
        //        template<typename... _Args>
        //        typed(_Args&&... args)
        //          : base_t(std::forward<_Args>(args)...)
        //        {}

        template<typename... Args> void init(Args&&... args)
        {
            Base::init([](base_jobtoken && token, wker & worker) noexcept {
                typed<F, ctx> ctx(std::move(token), worker);
                auto&         functor = ctx.functor();
                ctx.schedule_destructor();
                functor(std::move(ctx));
            });
            this->emplace(std::forward<Args>(args)...);
        }
    };

    template<typename Functor> struct subtask_functor : public Functor
    {
        using res_ty = decltype(std::declval<Functor>()(std::declval<ctx>()));
        template<typename... Args>
        subtask_functor(const jobtoken& parent, Args&&... args)
          : Functor(std::forward<Args>(args)...)
          , _parent(parent)
        {
            assert(_parent, "Task with empty parent");
        }
        const jobtoken _parent;
    };

    struct ctx : public jobtoken
    {
        ctx(jobtoken&& token, wker& worker)
          : jobtoken(token)
          , _worker(worker)
        {
            assert(token, "Context with empty token");
        }

        ctx(ctx&&)      = default;
        ctx(const ctx&) = delete;

        // Return a new reference to the job
        jobtoken token() const { return jobtoken(*this); }

        template<typename F, typename... Args> void subtask(Args&&... args)
        {
            _worker.schedule_job(_worker.allocate().template init<1, typed<subtask_functor<F>, jobtoken>>(
              *this, std::forward<Args>(args)...)[0]);
        }

        template<typename F> void subtask(F&& f) { return subtask<remove_reference_t<F>, F>(std::forward<F>(f)); }

        template<typename F, typename... Args> void schedule(Args&&... args)
        {
            _worker.schedule_job(
              _worker.allocate().template init<1, typed<F, jobtoken>>(*this, std::forward<Args>(args)...));
        }

        template<typename F> void schedule(F&& f) { return schedule<remove_reference_t<F>, F>(std::forward<F>(f)); }

        template<typename F, typename... Args> typed<subtask_functor<F>, jobtoken> subtask_with_token(Args&&... args)
        {
            auto tks = _worker.allocate().template init<2, typed<subtask_functor<F>, jobtoken>>(
              *this, std::forward<Args>(args)...);
            _worker.schedule_job(std::move(tks[0]));
            return std::move(tks[1]);
        }

        template<typename F> typed<subtask_functor<F>, jobtoken> subtask_with_token(F&& f)
        {
            return subtask_with_token<remove_reference_t<F>, F>(std::forward<F>(f));
        }

        template<typename F, typename... Args> typed<F, jobtoken> schedule_with_token(Args&&... args)
        {
            auto tks = _worker.allocate().template init<2, typed<F, jobtoken>>(*this, std::forward<Args>(args)...);
            _worker.schedule_job(std::move(tks[0]));
            return std::move(tks[1]);
        }

        template<typename F> typed<F, jobtoken> schedule_with_token(F&& f)
        {
            return schedule_with_token<remove_reference_t<F>, F>(std::forward<F>(f));
        }

        wker& _worker;
    };

    // template<typename Functor> using jobref = typename jobslot::template jobref<Functor>;
    const jobtoken& root_ref;
    jobtoken        root;

    // FIXME: emplace version
    template<typename F> static void start(unsigned ncpu, F&& f)
    {
        jobslot root;
        auto    tytks = root.template init<1, typed<F, jobtoken>>(std::forward<F>(f));

        assume(not root.empty(), "uninit");
        sys::run_pinned_worker_pool<wker>(ncpu, std::cref(tytks[0]));
    }

    wker(const jobtoken& root)
      : root(root)
      , root_ref(root)
      , _tid(std::this_thread::get_id())
    {
        assume(root_ref, "uninit");
    }

    void operator()(sys::worker_pool_t<wker> workers, unsigned worker_id)
    {
        order_siblings(workers, worker_id);

        if (worker_id == 0) {
            root.run(*this);
        } else {
            root.~jobtoken();
        }

        wait_untill([&]() { return root_ref.single_owner(); });
    }

  private:
    void order_siblings(sys::worker_pool_t<wker> workers, unsigned worker_id)
    {
        size_t nsiblings = workers.size() - 1;
        siblings         = make_unique<wker*[]>(nsiblings);

        // Sort siblings by distance to this worker, starting with the worker in the same pair as this one
        // (hyper-thread)
        for (unsigned i = 1; _nsiblings < nsiblings; ++i) {
            unsigned d = (i + (worker_id & 1)) >> 1;
            if (d == 0) continue;
            if (i & 1) {
                if (d <= worker_id) siblings[_nsiblings++] = workers[worker_id - d];
            } else {
                if (worker_id + d < workers.size()) siblings[_nsiblings++] = workers[worker_id + d];
            }
        }
    }

    jobslot& allocate() noexcept
    {
        jobslot* j = _pool.allocate();
        for (; unlikely(j == nullptr); j = _pool.allocate())
            yield();
        return *j;
    }

    void schedule_job(jobtoken j)
    {
        // Wait while the queue is full
        for (j = _jobs.push(std::move(j)); unlikely(j); j = _jobs.push(std::move(j)))
            yield();
    }

    template<typename P> forceinline_fun void wait_untill(P&& predicate)
    {
        while (not predicate())
            yield();
    }

    void yield(unsigned steal_retries = 1) noexcept
    {
        assert(check_owner(), "worker.get_job called from the wrong thread");

        jobtoken j = _jobs.pop();
        if (j) {
            j.run(*this);
            return;
        }

        j = steal(steal_retries);

        if (j) j.run(*this);
    }

    jobtoken steal(unsigned steal_retries = 1) const noexcept
    {

        for (unsigned i = 0; i < steal_retries; i++) {
            wker* busiest      = nullptr;
            int   busiest_njob = 0;
            for (unsigned i = 0; i < _nsiblings; i++) {
                wker* worker = siblings[i];
                int   njobs  = worker->_jobs.approximate_size();
                if (njobs > busiest_njob) {
                    busiest      = worker;
                    busiest_njob = njobs;
                }
            }
            if (busiest != nullptr)
                return busiest->_jobs.steal();
            else
                continue;
        }
        return jobtoken();
    }

    bool check_owner() const noexcept { return std::this_thread::get_id() == _tid; }

    std::thread::id          _tid;
    std::unique_ptr<wker*[]> siblings;
    unsigned                 _nsiblings{};

    scanning_pool_allocator<jobslot, pool_bits> _pool;
    work_stealing_queue<jobtoken, cap_bits>     _jobs = {};
};

} // namespace sched
} // namespace gatbl

#endif // SCHEDULER_HPP