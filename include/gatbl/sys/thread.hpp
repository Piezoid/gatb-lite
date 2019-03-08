#ifndef THREAD_HPP
#define THREAD_HPP

#include <pthread.h>
#include <stdlib.h> // For posix_memalign

#ifdef _GNU_SOURCE
#    include <sched.h> // for sched_getcpu
#endif

#include <vector>
#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "gatbl/sys/exceptions.hpp"

namespace gatbl { namespace sys {

/// Pin the current thread to a single cpu and wait for it's migration
inline void
pin_to_cpu(unsigned cpu)
{
#ifdef _GNU_SOURCE
    cpu_set_t cpuset;
    pthread_t thread;

    thread = pthread_self();

    CPU_ZERO(&cpuset);
    CPU_SET(cpu, &cpuset);
    sys::check_ret(pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset), "pthread_setaffinity_np failed");

    pthread_yield();
    for (unsigned i = 0; sys::check_ret(sched_getcpu(), "sched_getcpu") != cpu; ++i) {
        pthread_yield();
        if (i >= 1024) {
            std::cerr << "Failled to migrate to CPU " << cpu << " (currently on " << sched_getcpu()
                      << "). Is the system too busy ?\n";
            break;
        }
    }
#endif
}

template<typename Worker> using worker_pool_t = const std::vector<Worker*>&;

/// Construct and run a worker pool with threads pinned on CPUs.
///
/// When worker objects are allocated on the target cpu they receive a copy of arguments.
/// Once all workers are constructed, a synchronization point is reached and they are called with a vector of pointers
/// to the workers (alias worker_pool_t<Worker>) ordered by CPUs and their CPU number (0-based) Worker have the
/// responsibility of managing their running condition.
template<typename Worker, typename... Args>
inline void
run_pinned_worker_pool(unsigned ncpu, const Args&&... args)
{
    std::vector<Worker*>     workers{ncpu};
    std::vector<std::thread> threads;
    threads.reserve(ncpu);

    std::mutex              mut{};
    unsigned                started = 0;
    std::condition_variable cond{};

    auto make_worker = [&](unsigned thread_num) {
        // Thread are pinned and migrated to CPU before allocating the Worker data, enabling better NUMA locality under
        // default policy
        pin_to_cpu(thread_num);

        // FIXME: correct alignement for overaligned types is only respected for after cpp17
        // auto worker = make_unique<Worker>(Args(static_cast<const Args&>(args))...);
        Worker* worker;
        sys::check_ret(posix_memalign(reinterpret_cast<void**>(&worker), alignof(Worker), sizeof(Worker)),
                       "Failled allocating aligned memory");
        new (worker) Worker(Args(static_cast<const Args&>(args))...);

        {
            std::unique_lock<std::mutex> lock(mut);
            workers[thread_num] = worker;
            if (++started < ncpu)
                cond.wait(lock);
            else
                cond.notify_all();
            assume(started == ncpu, "Worker starting while others are not ready");
        }
        worker->operator()(worker_pool_t<Worker>(workers), thread_num);

        // FIXME: see previous fixme
        worker->~Worker();
        free(worker);
    };

    // Spawn threads for ncpu-1 workers
    for (unsigned i = 1; i < ncpu; i++)
        threads.emplace_back(make_worker, i);

    // Main thread worker
    make_worker(0);

    for (auto& thread : threads)
        thread.join();
}

} // namespace sys
} // namespace gatbl

#endif // THREAD_HPP
