
#include <cstddef>
#include <iostream>
#include <cstring>

#include "gatbl/kmer.hpp"
#include "gatbl/utils/string.hpp"

using namespace gatbl;
using namespace std;

int
main(int argc, char** argv)
{
    ksize_t k     = 15;
    using kmer_t  = uint64_t;
    using skmer_t = sized_kmer<kmer_t>;

    entropy_filter<kmer_t> filter(k, 0);

    constexpr size_t hist_size = filter.denominator() * 4;
    size_t           hist[hist_size];
    for (auto& x : hist)
        x = 0;

    kmer_t kmax = kmer_t(1) << (2 * k);
    for (kmer_t kmer = 0; kmer < kmax; kmer++) {
        size_t n = filter.numerator(kmer);
        assume(n < hist_size, "wtf");
        hist[n]++;
        // if (entropy < 1) std::cout << skmer_t{kmer, k} << ": " << filter.entropy(kmer) << std::endl;
    }

    for (auto& x : hist) {
        std::cout << x << " ";
    }
    std::cout << std::endl;
}
