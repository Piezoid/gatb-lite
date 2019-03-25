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
    ~bounded_size_wrapper() noexcept(noexcept(std::declval<value_type>().~value_type())) { delete value; }

  private:
    value_type* value;
};

} // namespace memory

namespace sched {

constexpr size_t hardware_destructive_interference_size = 64;

/** jobslot is a like a reference counted std::function
 *
 * The functor data is in inline storage, up to some number of cache lines.
 *
 * Users must derive the jobtoken (smart pointer reference counter) and typed (wrapper over jobtoken
 * derived classes for unforgetting the type of the functor).
 */
template<typename Fun, size_t cache_lines = 1> class jobslot;

template<typename... Args, size_t cache_lines>
class alignas(hardware_destructive_interference_size) jobslot<void(Args...), cache_lines>
{
  public:
    class jobtoken;

  protected:
    /// Type of the raw function. it receive the token from witch it was invoked.
    using func_ptr = void (*)(jobtoken, Args...) CPP17_NOEXCEPT;
    /// Type of the finish handler. Expected to be set by the function. Called when the reference count reaches 0.
    using finish_handler_ptr = void (*)(void*) CPP17_NOEXCEPT;

  private:
    static constexpr size_t data_size
      = cache_lines * hardware_destructive_interference_size - sizeof(std::atomic<size_t>) - sizeof(func_ptr);
    template<typename F> using data_wrapper = memory::bounded_size_wrapper<F, data_size>;

    alignas(hardware_destructive_interference_size) char _data[data_size];

    std::atomic<size_t> _state_idx{0};
    // _state_idx = 0 is empty job ready to be allocated in the pool
    // _state_idx = 1 is a job being finished
    // _state_idx = 1 + n is scheduled job with (n-1) alive references
    func_ptr _function = nullptr;

    /// Decrement the state counter
    void release() noexcept
    {
        // jobslot type invariants are checked here, now that the type is complete
        static_assert(sizeof(jobslot) == cache_lines * hardware_destructive_interference_size,
                      "job is the size of some cache line");
        static_assert(alignof(jobslot) % hardware_destructive_interference_size == 0,
                      "job is aligned on cache line boundaries");
        assert(reinterpret_cast<uintptr_t>(this) % hardware_destructive_interference_size == 0,
               "Unaligned job slot"); // RT version of the above

        auto state_idx = _state_idx.fetch_sub(1,
                                              // TSan doesn't handle well the conditional fence here
                                              __has_feature(thread_sanitizer) ? std::memory_order_acq_rel
                                                                              : std::memory_order_release);
        assume(state_idx > 1, "Double free job slot");
        if (state_idx == 2) // If it reaches 1 the job completion handler is executed
        {
            std::atomic_thread_fence(std::memory_order_acquire);
            // Publishing _state_idx == 1 allows to postpone the reuse of the pool's slot while running the
            // completion chain.
            if (_function != nullptr) {
                // Call the completion handler
                reinterpret_cast<finish_handler_ptr>(_function)(_data);
                if (DEBUG) _function = nullptr;
            }

            // Mark the job as finished so it can be reallocated.
            // Pretty sure that relaxed ordering is enough here: We don't mind if the allocator sees the
            // previously released value _state_idx=1 instead of 0. Sadly, TSan thinks otherwise...
            static constexpr std::memory_order mo
              = __has_feature(thread_sanitizer) ? std::memory_order_release : std::memory_order_relaxed;
            if (DEBUG) {
                assert(_state_idx.compare_exchange_strong(--state_idx, 0, mo, std::memory_order_relaxed),
                       "Race condition in job release()");
            } else {
                _state_idx.store(0, mo);
            }
        }
    }

  public:
    // The triaval constructor make an empty job in finised state
    jobslot() noexcept      = default;
    jobslot(const jobslot&) = delete;
    jobslot& operator=(const jobslot&) = delete;

#if DEBUG
    ~jobslot() noexcept { assert(_state_idx < 2, "non executed job destroyed"); }
#else
    ~jobslot() noexcept = default;
#endif

    bool empty(std::memory_order mo = std::memory_order_acquire) const noexcept { return _state_idx.load(mo) == 0; }

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
        /// Execute the job, and drop the reference as if we moved from this token
        void run(Args... args) noexcept
        {
            jobslot& s = slot();
            // Precondition : have a function defined and alive jobrefs
            assert(s._state_idx.load(std::memory_order_acquire) >= 2, "Trying to run a disowned job");
            assert(*this, "Trying to run an empty job");
            s._function(std::move(*this), std::forward<Args>(args)...);
        }

        /// A very basic initializer for the job for raw function pointers
        /// NB: the _function must call onfinish before returning !
        void init(func_ptr f) { slot()._function = f; }

        /// Set the finish/join handler of the job. It is executed when the last reference to the job is destructed.
        /// This must be called before the job's function returns
        /// If the function is null, no action is taken.
        void onfinish(finish_handler_ptr finish = nullptr) { slot()._function = reinterpret_cast<func_ptr>(finish); }

      public:
        /// An empty jobref not referencing any jobconst
        jobtoken() noexcept = default;

        /// Duplicate ownership of the job
        explicit jobtoken(const jobtoken& other) noexcept
          : _j(other._j)
        {
            if (operator bool()) {
                size_t state_idx = _j->_state_idx.fetch_add(1, std::memory_order_relaxed);
                assume(state_idx >= 2, "acquire() on unititialized job");
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
            if (operator bool()) _j->release();
            _j       = other._j;
            other._j = nullptr;
            return *this;
        }

        /// Release the ownership of the job
        forceinline_fun ~jobtoken() noexcept
        {
            if (operator bool()) {
                _j->release();
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
            assert(operator bool(), "Asking owner count of a null token");
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

    /// Typed interface over a job token: allows to access and mutate the slot content as a Functor instance
    template<typename F, typename Base = jobtoken> class typed : public Base
    {
        using data_t = data_wrapper<F>;
        friend class jobslot;

      protected:
        /// Downcasting from an untyped token: as dangerous as reinterpret_cast !
        template<typename... _Args>
        typed(_Args&&... args)
          : Base(std::forward<_Args>(args)...)
        {}

        template<typename... ConstrArgs> void emplace(ConstrArgs&&... constr_args)
        {
            new (Base::slot()._data) data_t(std::forward<ConstrArgs>(constr_args)...);
        }

        /// Set deletion handler (should only be used once the job is started (e.g. at epilogue of the job function)
        /// Additionally a join handler can be specified as a member function, called before destruction
        template<void (F::*finish_member_ptr)() = nullptr> void onfinish()
        {
            finish_handler_ptr finish = nullptr;
            if (finish_member_ptr != nullptr || not std::is_trivially_destructible<data_t>::value)
                finish = [](void* ptr) noexcept
                {
                    auto& data = *reinterpret_cast<data_t*>(ptr);
                    if (finish_member_ptr != nullptr) { (data.get().*finish_member_ptr)(); }
                    if (not std::is_trivially_destructible<data_t>::value) { data.~data_t(); }
                };
            Base::onfinish(finish);
        }

        /// An example of job initialization that can be achieved with the above mutators
        template<typename... ConstrArgs> void init(ConstrArgs&&... args)
        {
            Base::init([](jobtoken token, Args... args) noexcept {
                // In the function, the type of the functor is forgoten, so get it back by constructing a typed token of
                // the right type
                typed tytk(std::move(token));
                auto& functor = tytk.functor();
                tytk.onfinish();
                functor(std::move(tytk));
            });
            this->emplace(std::forward<ConstrArgs>(args)...);
        }

      public:
        static constexpr bool have_indirection = data_wrapper<F>::have_indirection;

        F& functor() { return reinterpret_cast<data_wrapper<F>*>(Base::slot()._data)->get(); }

        const F& functor() const { return reinterpret_cast<data_wrapper<F>*>(Base::slot()._data)->get(); }

        using Base::Base;
    };

    /// Initialize the jobslot and return a token
    template<typename Token, typename... ConstrArgs>
    auto init(ConstrArgs&&... constr_args) -> decltype(Token(jobtoken()))
    {
        // State transition to an occupied jobslot
        assert(empty(), "Allocated a non freed job");
        _state_idx.store(2, std::memory_order_relaxed);

        Token token{};
        token.set_noacquire(this);
        token.init(std::forward<ConstrArgs>(constr_args)...);
        return token;
    }

    /// Initialize the jobslot and return a pair of token
    template<typename Token, typename... ConstrArgs>
    auto init_pair(ConstrArgs&&... constr_args) -> std::pair<decltype(Token(jobtoken())), Token>
    {
        // State transition to an occupied jobslot
        assert(empty(), "Allocated a non freed job");
        _state_idx.store(3, std::memory_order_relaxed);

        std::pair<Token, Token> token{};
        token.first.set_noacquire(this);
        token.first.init(std::forward<ConstrArgs>(constr_args)...);
        token.second.set_noacquire(this);
        return token;
    }
};

/** Concurrent SPMC stack/queue for work stealing
 *
 * The producer-owner can push and pop in FIFO order from one end while other thread might steal elements in LIFO order
 * on the other end. _top and _bottom atomic indexes are read and write positions of a ring buffer exploiting the
 * modular arithmetic of the unsigned indice type.
 * @see https://www.snellman.net/blog/archive/2016-12-13-ring-buffers/
 */
template<typename T, size_t cap_bits> struct work_stealing_queue
{
    using idx_t                     = size_t;
    static constexpr idx_t capacity = 1ULL << cap_bits;
    static constexpr idx_t mask     = capacity - 1;

    /** Try to push a item (reserved to the single-producer/owner)
     * In case of success (queue is not full) the ownership of the item is passed and the reference is left in the
     * moved-from state.
     */
    void push(T& job) noexcept
    {
        assert(job, "pushing empty job");
        // Bottom: relaxed load, since _bottom is only modified by the worker owning the queue
        idx_t b = _bottom.load(std::memory_order_relaxed); // Will be incremented

        // Top: relaxed load, since _top is monotonically increasing we only risk false
        // positives in the "is full" test bellow
        idx_t t = _top.load(std::memory_order_relaxed);

        if (likely(b - t < capacity)) { // not full
            new (&at(b++)) T(std::move(job));
            _bottom.store(b, std::memory_order_release);
        }
    }

    /// Pops in FIFO order (reserved to the single-producer/owner)
    T pop() noexcept
    {
        // First we decrement _bottom: load is relaxed since we are the only thread to update it,
        // but publishes the updated value for stealers
        idx_t bottom = _bottom.fetch_sub(1, std::memory_order_release);
        idx_t top    = _top.load(std::memory_order_acquire);

        T item{};

        if (likely(bottom - top > 1)) { // More than one job left
            item = std::move(at(bottom - 1));
            assume(item, "poping invalid item");
            return item;
        }

        if (bottom - top == 1) {
            // Last job in the queue: it could also be extracted concurrently from the other end by steal()
            // So we also do a kind of steal() here, increasing top with CMPXCHG
            if (_top.compare_exchange_weak(top, top + 1, std::memory_order_acq_rel, std::memory_order_relaxed)) {
                // we were alone incrementing top
                item = std::move(at(bottom - 1)); // Or at(top)
                assume(item, "poping invalid item");
            }
        }

        // restore bottom to the original value (we either self-stole, or found and empty queue)
        _bottom.store(bottom, std::memory_order_release);
        assert(_bottom == _top, "Unexpected top idx change");
        return item;
    }

    /// Steal an element in LIFO order. Can be called from any thread.
    T steal() noexcept
    {
        // Relaxed load since the value will be checked with a CAS
        idx_t t = _top.load(std::memory_order_relaxed); // Will be incremented
        idx_t b = _bottom.load(std::memory_order_acquire);

        // The queue not empty if bottom > top. Howerver, we can't test that because we want to tolerate wrapparound of
        // the counter.
        // empty => bottom == top, unless bottom is temporarly decremented by the queue owner inspecting it
        // with pop() in wich case t == b + 1.
        if (t != b && t != b + 1) {
            // FIXME: nasty lookahead. We need to move the object before releasing the incremented top, but without the
            // side effects of moving in case of backtracking.
            storage_t item_raw = _arr[t & mask];
            if (_top.compare_exchange_weak(t, t + 1, std::memory_order_acq_rel, std::memory_order_relaxed)) {
                T& item = reinterpret_cast<T&>(item_raw);
                assume(item, "Stealed an empty item");
                return std::move(item);
            } else {
                return T{};
            }
        }

        return T{};
    }

    int approximate_size() const noexcept
    { // TODO: check if its a win to avoid a memory fence here (might loose oportunities for stealing)
        return int(_bottom.load(std::memory_order_relaxed)) - int(_top.load(std::memory_order_relaxed));
    }

#if DEBUG
    ~work_stealing_queue() { assert(_bottom == _top, "Jobs remaining in the queue upon destruction"); }
#endif

  private:
    T& at(idx_t i) { return reinterpret_cast<T&>(_arr[i & mask]); }

    using storage_t = typename std::aligned_storage<sizeof(T), alignof(T)>::type;

    storage_t _arr[capacity];
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

  private:
    static constexpr size_t _mask = size - 1;
    value_type              _pool[size];
    size_t                  _idx = 0;
};

/**
 * A worker, with a start() static method for initiating a worker pool executing subtasks spawn from a single root task.
 *
 * Work in progress, API might change a lot.
 */
template<size_t cap_bits = 6, size_t pool_bits = cap_bits> class worker
{
    using base_jobslot                                   = jobslot<void(worker&)>;
    using base_jobtoken                                  = typename base_jobslot::jobtoken;
    template<typename F, typename Base> using base_typed = typename base_jobslot::template typed<F, Base>;

  public:
    struct ctx;
    struct jobtoken : public base_jobtoken
    {
        using Base = typename base_jobslot::jobtoken;
        using base_jobtoken::base_jobtoken;

      protected:
        template<typename Tk>
        jobtoken(Tk&& from)
          : base_jobtoken(std::forward<Tk>(from)){};

        friend class worker;
        friend base_jobslot;
    };

    template<typename F, typename Base = jobtoken> struct typed : public base_typed<F, Base>
    {
        using base_t = base_typed<F, Base>;
        using base_t::base_t;

        // Return a new reference to the job
        typed<F, jobtoken> token() const { return typed<F, jobtoken>(*this); }

        // Return an upcasted version of the same reference
        typed<F, jobtoken> astoken() const { return typed<F, jobtoken>(std::move(*this)); }

      protected:
        friend class worker;
        friend base_jobslot;
        //        template<typename... _Args>
        //        typed(_Args&&... args)
        //          : base_t(std::forward<_Args>(args)...)
        //        {}

        template<typename... Args> void init(Args&&... args)
        {
            Base::init([](base_jobtoken token, worker & worker) noexcept {
                assume(token, "Job called with from an empty token");
                typed<F, ctx> ctx(std::move(token), worker);
                auto&         functor = ctx.functor();
                ctx.onfinish();
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
            assume(_parent, "Task with empty parent");
        }
#ifdef DEBUG
        // In case the job slot is overwritten
        ~subtask_functor() { assume(_parent, "Task with empty parent"); }
#endif
        const jobtoken _parent;
    };

    struct ctx : public jobtoken
    {
        ctx(jobtoken&& token, worker& worker)
          : jobtoken(std::move(token))
          , _worker(worker)
        {
            assert(jobtoken::operator bool(), "Context with empty token");
        }

        ctx(ctx&&)      = default;
        ctx(const ctx&) = delete;

        // Return a new reference to the job
        jobtoken token() const { return jobtoken(*this); }

        template<typename F, typename... Args> void subtask(Args&&... args)
        {
            _worker.schedule_job(_worker.allocate().template init<typed<subtask_functor<F>, jobtoken>>(
              *this, std::forward<Args>(args)...));
        }

        template<typename F> void subtask(F&& f) { return subtask<remove_reference_t<F>, F>(std::forward<F>(f)); }

        template<typename F, typename... Args> void schedule(Args&&... args)
        {
            _worker.schedule_job(
              _worker.allocate().template init<typed<F, jobtoken>>(*this, std::forward<Args>(args)...));
        }

        template<typename F> void schedule(F&& f) { return schedule<remove_reference_t<F>, F>(std::forward<F>(f)); }

        template<typename F, typename... Args> typed<subtask_functor<F>, jobtoken> subtask_with_token(Args&&... args)
        {
            auto tks = _worker.allocate().template init_pair<typed<subtask_functor<F>, jobtoken>>(
              *this, std::forward<Args>(args)...);
            _worker.schedule_job(std::move(tks.first));
            return std::move(tks.second);
        }

        template<typename F> typed<subtask_functor<F>, jobtoken> subtask_with_token(F&& f)
        {
            return subtask_with_token<remove_reference_t<F>, F>(std::forward<F>(f));
        }

        template<typename F, typename... Args> typed<F, jobtoken> schedule_with_token(Args&&... args)
        {
            auto tks = _worker.allocate().template init_pair<typed<F, jobtoken>>(*this, std::forward<Args>(args)...);
            _worker.schedule_job(std::move(tks.first));
            return std::move(tks.second);
        }

        template<typename F> typed<F, jobtoken> schedule_with_token(F&& f)
        {
            return schedule_with_token<remove_reference_t<F>, F>(std::forward<F>(f));
        }

        template<typename P> forceinline_fun void wait_untill(P&& predicate)
        {
            _worker.wait_untill(std::forward<P>(predicate));
        }

      private:
        worker& _worker;
    };

    jobtoken root; // Token to the root job, only live during the time between construction and execution of the worker.
    const jobtoken& root_ref; // Reference to the token on the stack of the thread starting the worker pool, live for
                              // the whole work pool lifetime

    // FIXME: emplace version
    template<typename F> static void start(unsigned ncpu, F&& f)
    {

        base_jobslot root; // Slot holding the main task
        // This token, held during the lifetime of the work pool, act as a witnedd for the termination of all the tasks
        const auto token = root.template init<typed<F, jobtoken>>(std::forward<F>(f));
        sys::run_pinned_worker_pool<worker>(ncpu, std::cref(token));
    }

    worker(const jobtoken& root)
      : root(root)
      , root_ref(root)
      , _tid(std::this_thread::get_id())
    {
        assert(root_ref, "Root task is null");
    }

    void operator()(sys::worker_pool_t<worker> workers, unsigned worker_id)
    {
        order_siblings(workers, worker_id);

        // Each worker have a reference increment to the root job untill now
        if (worker_id == 0) {
            // Consume the token, but might generate new reference in sub tasks
            root.run(*this);
        } else {
            // Delete our reference since worker0 handle it's execution
            root = jobtoken();
        }

        // Run jobs until there is only one token reference to the jobslot on start()'s stack frame
        wait_untill([&]() { return root_ref.single_owner(); });
    }

  private:
    void order_siblings(sys::worker_pool_t<worker> workers, unsigned worker_id)
    {
        size_t nsiblings = workers.size() - 1;
        siblings         = make_unique<worker*[]>(nsiblings);

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

    base_jobslot& allocate() noexcept
    {
        base_jobslot* j = nullptr;
        wait_untill([&]() {
            j = _pool.allocate();
            return j != nullptr;
        });
        assert(j->empty(), "Allocated a non empty job");
        return *j;
    }

    void forceinline_fun hot_fun schedule_job(base_jobtoken j)
    {
        wait_untill([&]() {
            _jobs.push(j);
            return not(j);
        });
    }

    template<typename P> forceinline_fun void wait_untill(P&& predicate)
    {
        assert(check_owner(), "worker.get_job called from the wrong thread");

        if (predicate())
            return;
        else {
            std::pair<jobtoken, size_t> res;
            res.second = 0;
            do {
                res = get_job(res.second);
                if (res.second == 0) res.first.run(*this);
            } while (not predicate());
        }
    }

    // In order to avoid code bloat, the code for getting a job is behind a call
    // However, the job is run on the parent stack frame.
    noinline_fun hot_fun std::pair<jobtoken, size_t> get_job(size_t state = 0)
    {
        // If this is the first try, we pop from our queue
        if (likely(state == 0)) {
            jobtoken j = _jobs.pop();
            if (j)
                return {std::move(j), 0};
            else
                ++state; // Goes on to stealling
        } else
            assert(_jobs.approximate_size() == 0, "Asked to steal while we still have local jobs");

        // If this is the second *steal* try, we yield to the OS before retry
        if (state > 1) sched_yield();

        // If poping from our queue failed, we steal from other workers (in locality order)
        for (unsigned i = 0; i < _nsiblings; i++) {
            jobtoken j = siblings[i]->_jobs.steal();
            if (j) return {std::move(j), 0};
        }

        // Everything failed, return with state > 1
        return {jobtoken(), ++state};
    }

    bool check_owner() const noexcept { return std::this_thread::get_id() == _tid; }

    std::thread::id            _tid;
    std::unique_ptr<worker*[]> siblings;
    unsigned                   _nsiblings{};

    scanning_pool_allocator<base_jobslot, pool_bits> _pool{};
    work_stealing_queue<base_jobtoken, cap_bits>     _jobs{};
};

} // namespace sched
} // namespace gatbl

#endif // SCHEDULER_HPP
