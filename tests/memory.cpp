#include "gatbl/sys/memory.hpp"
#include "gatbl/ds/rings.hpp" //FIXME: add tests

#include <catch2/catch.hpp>

using namespace gatbl;

struct InstanceCounter
{
    InstanceCounter(uint8_t byte)
      : byte(byte)
    {
        count++;
    }
    InstanceCounter()
      : byte(0xaa)
    {
        count++;
    }
    ~InstanceCounter() { count--; }
    static int count;
    uint8_t    byte;
};

int InstanceCounter::count = 0;

struct alignas(512) InstanceCounterOverAligned : InstanceCounter
{
    using InstanceCounter::InstanceCounter;
};

template<typename F> void inline check_instance_count(F&& f)
{
    REQUIRE(InstanceCounter::count == 0);
    f();
    WHEN("unique_ptr goes out of scope")
    {
        THEN("Instance count returns to 0") { CHECK(InstanceCounter::count == 0); }
    }
}

SCENARIO("InstanceCounter sanity check", "[memory]")
{
    check_instance_count([]() {
        InstanceCounter inst{};
        REQUIRE(InstanceCounter::count == 1);
    });
}

template<typename T>
void
check_ptr(const T* ptr)
{
    THEN("The pointer is not null") { REQUIRE(ptr != nullptr); }

    THEN("Memory is correclty aligned to " << alignof(T) << "bytes")
    {
        CHECK((reinterpret_cast<uintptr_t>(ptr) & uintptr_t(alignof(T) - 1)) == 0);
    }

    THEN("The constructor was called") { CHECK(ptr->byte == 0xaa); }
}

template<typename T>
void
do_aligned_alloc_test()
{
    WHEN("Allocating a single object")
    {
        check_instance_count([]() { check_ptr(make_unique_aligned<T>(uint8_t{0xaa}).get()); });
    }

    WHEN("Allocating an array of 16 object")
    {
        check_instance_count([]() { check_ptr(make_unique_aligned<T[]>(16).get()); });
    }
}

SCENARIO("make_unique_aligned allocation", "[memory]")
{
    GIVEN("a non triavlly constructible/destructible type") { do_aligned_alloc_test<InstanceCounter>(); }

    GIVEN("a non triavlly constructible/destructible over-aligned type")
    {
        do_aligned_alloc_test<InstanceCounterOverAligned>();
    }
}
