#include <vector>
#include <catch2/catch.hpp>

#include "gatbl/sys/file.hpp"
#include "gatbl/ds/int_vector.hpp"

#include "test_utils.hpp"

#include <iostream> // FIXME
#include <fstream>

using namespace gatbl;

template<typename LHS, typename RHS>
void
check_equal(const LHS& lhs, const RHS& rhs)
{
    REQUIRE(lhs.size() <= rhs.size());
    THEN("It compare equal to the original vector") { CHECK(std::equal(lhs.begin(), lhs.end(), rhs.begin())); }
}

template<typename V>
void
check_size(const V& vec, size_t size)
{
    THEN("it's size is" << size)
    {
        CHECK(std::distance(vec.begin(), vec.end()) == std::ptrdiff_t(size));
        CHECK(vec.size() == size);
    }
}

template<typename Vec = int_vector<>>
void
check_int_vector(size_t sz = 4096, typename Vec::width_t width = 8 + 4, typename Vec::value_type val_init = 0xaa)
{
    using type = typename Vec::value_type;

    GIVEN("A vector type " << type_name<Vec>() << ", width=" << size_t(width) << ", size=" << sz)
    {

        GIVEN("A vector of integers of width initialized to " << val_init)
        {
            Vec vec(sz, width, val_init);
            THEN("All element compare equal to " << val_init)
            CHECK(std::all_of(vec.begin(), vec.end(), [&](type x) { return x == val_init; }));
        }

        GIVEN("A vector random integers")
        {
            auto stdvec = make_randint_vector<type>(sz, width);
            Vec  vec(stdvec);
            check_equal(vec, stdvec);

            size_t new_size = (2 * sz) / 3;
            WHEN("sized-down to " << new_size)
            {
                const size_t original_size_in_bytes = vec.size_in_bytes();
                vec.resize(new_size);

                check_size(vec, new_size);
                check_equal(vec, stdvec);

                THEN("Its size_in_bytes doesn't change") { CHECK(vec.size_in_bytes() == original_size_in_bytes); }

                AND_WHEN("shrinked to fit")
                {
                    vec.shrink_to_fit();
                    THEN("It's size_in_bytes is reduced") { CHECK(vec.size_in_bytes() < original_size_in_bytes); }
                    check_size(vec, new_size);
                    check_equal(vec, stdvec);
                }
            }

            new_size = (4 * sz) / 3;
            WHEN("sized-up to " << new_size)
            {
                vec.resize(new_size);
                check_size(vec, new_size);
                check_equal(stdvec, vec);
            }

            WHEN("assigned")
            {
                Vec copy;
                copy = vec;
                check_size(copy, sz);
                check_equal(copy, stdvec);
            }

            WHEN("copyed")
            {
                Vec copy(vec);
                check_size(copy, sz);
                check_equal(copy, stdvec);
            }

            WHEN("move-assigned")
            {
                Vec assigned;
                assigned = std::move(vec);
                check_size(assigned, sz);
                check_equal(assigned, stdvec);
            }

            WHEN("moved from")
            {
                Vec moved(std::move(vec));
                check_size(moved, sz);
                check_equal(moved, stdvec);
            }

            WHEN("sorted")
            {
                std::sort(stdvec.begin(), stdvec.end());
                std::sort(vec.begin(), vec.end());
                check_equal(vec, stdvec);
            }
        }
    }
}

SCENARIO("int_vector operations", "[int_vector]")
{
    check_int_vector<bit_vector>(4096, {}, true);
    check_int_vector<twobit_dna_vector>(2048, {}, nuc_t(2));
    check_int_vector<int_vector<packedint::specification<uint16_t, uint8_t, uint32_t>>>(512, 9, uint16_t(0x154));
    check_int_vector<int_vector<>>(256, 32, 0xdeadbeef);
    check_int_vector<int_vector<>>(1369, 11, 0x555);
    check_int_vector<int_vector<>>(128, 57, 0xfedcbadeadbeef);
}

template<typename Vec = int_vector<>>
void
check_int_vector_serialization(size_t sz = 4096, typename Vec::width_t width = 8 + 4)
{
    using std::ios;
    using type = typename Vec::value_type;

    GIVEN("A vector type " << type_name<Vec>() << ", width=" << size_t(width) << ", size=" << sz)
    {

        const char* filename = "example.bin";
        auto        stdvec   = make_randint_vector<type>(sz, width);

        // Serialize using std::fstream
        {
            std::ofstream file;
            file.open(filename, ios::out | ios::trunc | ios::binary);
            Vec vec(stdvec);
            write(file, vec);
        }

        {
            std::ifstream file;
            file.open(filename, ios::in | ios::binary);
            Vec vec{};
            read(file, vec);
            check_equal(stdvec, vec);
        }

        // Same but using file_descriptor
        {
            file_descriptor fd(filename, O_RDWR, O_CREAT | O_TRUNC);
            Vec             vec(stdvec);
            write(fd, vec);
        }

        {
            file_descriptor fd(filename);
            Vec             vec{};
            read(fd, vec);
            check_equal(stdvec, vec);
        }
    }
}

SCENARIO("int_vector serialization", "[int_vector][serialization]")
{
    check_int_vector_serialization<bit_vector>(4096, {});
    check_int_vector_serialization<twobit_dna_vector>(2048, {});
    check_int_vector_serialization<int_vector<packedint::specification<uint16_t, uint8_t, uint32_t>>>(512, 9);
    check_int_vector_serialization<int_vector<>>(256, 32);
    check_int_vector_serialization<int_vector<>>(128, 57);
}
