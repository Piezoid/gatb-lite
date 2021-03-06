#include <catch2/catch.hpp>

#include "gatbl/sys/file.hpp"
#include "gatbl/fastx.hpp"
#include "test_utils.hpp"

using namespace gatbl;

template<typename R>
inline void
count_test_records(size_t& count, const R& range, const size_t correct_length = 2)
{
    RANGES_FOR(auto rec, range)
    {
        CHECK(rec.sequence().size() == correct_length);
        count++;
    }
}

template<typename record_t, typename R>
inline void
sequence_iter_test(const R& data, const size_t correct_seq_num = 4, const size_t correct_length = 2)
{
    GIVEN("A char range of type " << type_name<R>() << " and a record type " << type_name<record_t>())
    {

        using seqit_t     = sequence_iterator<record_t>;
        using seq_range_t = sequence_range<record_t>;

        GIVEN("A full range sequence iterator: " << type_name<seq_range_t>())
        {
            seq_range_t r = data;
            THEN("Correct number of sequences of correct length are found")
            {
                size_t n = 0;
                count_test_records(n, r, correct_length);
                CHECK(n == correct_seq_num);
            }
        }

        for (size_t i = 0; i < size(data); ++i) {
            GIVEN("A range split at the " << i << "th char into two" << type_name<seq_range_t>())
            {
                seqit_t     mid_point(begin(data) + i, end(data));
                seq_range_t part1{begin(data), mid_point};
                seq_range_t part2{mid_point, end(data)};

                THEN("Correct number of sequences of correct length are found")
                {
                    size_t n = 0;
                    count_test_records(n, part1, correct_length);
                    count_test_records(n, part2, correct_length);
                    CHECK(n == correct_seq_num);
                }
            }
        }
    }
}

template<typename record_t>
inline void
mmap_sequence_iter_test(const char* filename, const size_t correct_seq_num = 4, const size_t correct_length = 2)
{
    GIVEN("A mmap of this file")
    {
        file_descriptor fd(filename);
        auto            data = fd.mmap<const char>();

        sequence_iter_test<record_t>(data, correct_seq_num, correct_length);
    }
}

template<typename record_t>
inline void
file_sequence_iter_test(const char* filename, const size_t correct_seq_num = 4, const size_t correct_length = 2)
{
    GIVEN("An input test file" << filename)
    {
        mmap_sequence_iter_test<record_t>(filename, correct_seq_num, correct_length);

        // FIXME: add other way of reading files here
    }
}

SCENARIO("FASTA file reading", "[IO][FASTX]") { file_sequence_iter_test<fasta_record<>>("data/test.fa"); }

SCENARIO("FASTQ file reading", "[IO][FASTX]") { file_sequence_iter_test<fastq_record<>>("data/test.fq"); }
