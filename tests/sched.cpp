

#include <algorithm>
#include <memory>
#include <vector>

#include <chrono>
#include <thread>
#include <pthread.h>
#include <cpuid.h>
#include <immintrin.h> // for _mm_pause()

#include "signal.h" //FIXME: debug

#include <stdlib.h> // For posix_memalign

#include <sched.h>

#include <iostream>

#include <execinfo.h>
void
print_backtrace()
{
    void* array[24];
    int   size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 24);
    backtrace_symbols_fd(array, size, 2);
}

#include "catch2/catch.hpp"
#include "gatbl/scheduler.hpp"
#include "gatbl/sys/exceptions.hpp"
#include "test_utils.hpp"

// using namespace gatbl::mem;
std::mutex global_out;
template<typename F>
inline void
with_mut(F&& f)
{
    std::lock_guard<std::mutex> lock(global_out);
    f();
}

using namespace gatbl;
using namespace gatbl::sched;

SCENARIO("start wp")
{

    using worker_t = worker<2, 5>;

    std::atomic<unsigned> counter{0};
    const unsigned        target_counts = 1024 * 1024; // std::numeric_limits<unsigned>::max() / 1024;

    // for (unsigned j = 0; j < 1024; j++)
    worker_t::start(8, [&](worker_t::ctx ctx) {
        std::cerr << "Started" << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        assert(ctx.owners_count() == 2, "%lu", ctx.owners_count());
        for (unsigned i = 0; i < target_counts; i++)
            ctx.subtask([&](worker_t::ctx) {
                // std::this_thread::sleep_for(std::chrono::milliseconds(100));
                counter++;
            });

        std::cerr << "Done" << std::endl;

        // ctx.template subtask<testfunctor>();
        // auto tk = ;
    });
    std::cerr << "All Done" << std::endl;

    CHECK(counter == target_counts);
}

/*
SCENARIO("job creation and destruction", "[scheduler]") {

    //static_assert(std::is_pod_v<job<>>, "job type should be a POD");

    // Make a job owning a shared pointer while we only keep a weak pointer to the same object.
    auto shared_j1 = std::make_shared<int>(0);
    auto shared_j2 = std::make_shared<int>(0);
    uint64_t test = 1;

    bool ran_j1 = false, ran_j2 = false;
    job<> j1 = [shared_j1](job<>&){};
    // The second job is a child of j1. The liftime of j1 is extended by j2
    job<> j2 = { [shared_j2](job<>&){ }, &j1 };

    // Disown the pointers after taking a weak reference
    auto weak_j1 = std::weak_ptr<int>(shared_j1); shared_j1.reset();
    auto weak_j2 = std::weak_ptr<int>(shared_j2); shared_j2.reset();

    REQUIRE(!weak_j1.expired()); REQUIRE(!weak_j2.expired());

    // Run the first job
    j1();

    // j1 is still alive since j2 is
    REQUIRE(!weak_j1.expired()); REQUIRE(!weak_j2.expired());

    j2();

    // Now they are both destroyed
    REQUIRE(weak_j1.expired()); REQUIRE(weak_j2.expired());
}
*/

/*

using worker_t = worker<>;

using worker_pool_t = worker_pool<worker_t>;
using job_t = typename worker_t::job_t;

#include<atomic>
#include <chrono>


template<typename T=uint64_t>
std::unique_ptr<uint64_t[]> binomials(uint64_t n0) {
    uint64_t g, r = 1;
    uint64_t n = n0;

    std::unique_ptr<uint64_t[]> arr(new uint64_t[n0 + 1]);

    arr[0] = r;
    for (uint64_t d = 1; d <= n0; d++) {
        r *= n--;
        r /= d;

        arr[d] = r;
    }

    return arr;
}


using atomic_array = std::unique_ptr<std::atomic<uint64_t>[]>;

//struct bd_ {

//    void operator() (job_t& j, worker_t& w, ctx&) const noexcept {

//        if(depth < n) {
//            w.schedule(bd_{arr, depth+1, n, k});
//            w.schedule(bd_{arr, depth+1, n, k+1}, j);
//            //std::this_thread::yield();
//            //w.wait_untill([&](){ return left.finished() && right.finished(); });
//        } else {
//            arr[k].fetch_add(1, std::memory_order_relaxed);
//        }
//    }

//    atomic_array& arr;
//    uint64_t depth, n, k;
//};

#include <functional>

using namespace std::chrono_literals;


std::function<void(job_t& j, worker_t& w)> global_fun;

struct bd {


    bd(atomic_array& arr, uint64_t depth, uint64_t n, uint64_t k) :
        arr(arr), depth(depth), n(n), k(k) {}

    bd(bd&&) noexcept = default;
    //bd(bd&) noexcept = default;

    template<typename Tok>
    void operator() (Tok tok) const noexcept {
//        for(unsigned i = 0 ; i < this->depth ; i++) std::cout << " ";
//        std::cout << depth << " " << k << ": started on " << std::this_thread::get_id() << std::endl;

        for(auto d =  depth ; d < n ; d++) {
            auto jt = job_token<job_t>(tok._job.get());
            tok._worker.template subschedule<bd>(std::move(jt), arr, d+1, n, k+1);
            //global_fun = std::function<void(job_t& j, worker_t& w, ctx&)>(bd{arr, d+1, n, k+1, j});
        }


//        std::this_thread::sleep_for(1ms);
//        for(unsigned i = 0 ; i < this->depth ; i++) std::cout << " ";
//        std::cout << depth << " " << k << ": done"  << std::endl;
    }

    ~bd() {
//        for(unsigned i = 0 ; i < this->depth ; i++) std::cout << " ";
//        std::cout << depth << " " << k << ": deleted on " << std::this_thread::get_id() << std::endl;
        arr[k].fetch_add(1, std::memory_order_relaxed);
    }

    atomic_array& arr;
    uint64_t depth, n, k;
    char x[64];
};


SCENARIO("worker pool", "[scheduler]") {

    const uint64_t n = 20;

    std::unique_ptr<std::atomic<uint64_t>[]> arr(new std::atomic<uint64_t>[n + 1]);
    for(uint64_t i = 0 ; i <= n ; i++) arr[i] = 0;

    auto ref = binomials(n);

    worker_pool_t wp(8);

    worker_t& me = wp.get_worker();




    auto start = std::chrono::steady_clock::now();
    me.schedule_wait<bd>(arr, 0, n, 0);
    auto end = std::chrono::steady_clock::now();
    wp([](){ return true;});


    //assume(j->finished());


    uint64_t sum = 0;
    for(uint64_t i = 0 ; i <= n ; i++) {
        REQUIRE(arr[i] == ref[i]);
        sum += arr[i];
    }
    std::cout << "Number of job executed:" << sum << ", time per job:"
              << std::chrono::duration <double, std::nano> ((end - start) / sum ).count() << " ns\n";
}

//SCENARIO("worker pool", "[scheduler]") {



//    std::atomic<unsigned> counter = 0;
//    worker_pool_t wp(32);

//    worker_t& me = wp.get_worker();

//    auto& j = me.schedule([&](job_t& j, worker_t& w, ctx&){
//        const auto fun = [&](job_t& j, worker_t& w, ctx&){
//            counter += 1;
//            std::this_thread::yield();
//        };
//        for(unsigned i = 0 ; i < 4096 ; i++) {
//            std::this_thread::yield();
//            w.schedule(fun, j);
//        }
//    });

//    wp([&]() { return j.finished(); });
//    //me.wait_untill([&]() { return j.finished(); });
//    REQUIRE(counter == 4096);


//    binomials(5);

//    //job_pool_allocator<job_t> pool;

*/
//}
