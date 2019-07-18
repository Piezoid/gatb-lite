
#include <cstddef>
#include <iostream>
#include <cstring>

#include "gatbl/fastx2kmer.hpp"

using namespace gatbl;
using namespace std;

int
main(int argc, char** argv)
{
    using kmer_t    = uint64_t;
    const ksize_t k = 16; // Kmer size. ]0, 32] for kmer_t = uint64_t

    size_t id_seq = 0;

    // The high level class instance that handle parsing fast[aq], filling the window, and emptying it on new sequence
    // or invalid nucleotides.
    auto f2kmers = make_sequence2kmers( //
      kmer_window<kmer_t>(k),           // The class handling the shifting of nucleotides into kmers
      // List of callback for different events. If ommited, they default to no operation.
      // The argument to CBs are the return value of make_sequence2kmers, that's the only way to have them in scope
      // (needs C++14 for auto lambda arguments)
      [&](auto& f2kmers) { // Callback for new kmer
          sized_kmer<kmer_t> kmer = f2kmers.get_window().forward();
          std::cout << kmer // sized_kmer is printable since it know it real size
                    << " " << f2kmers.get_pos() << std::endl;
      },
      [&](auto& f2kmers) { // New sequence in fast[aq]
          std::cout << "sequence " << id_seq++ << " start:" << f2kmers.get_next_pos() << std::endl;
      },
      [&](auto& f2kmers) { // Ends
          std::cout << "end run " << f2kmers.get_pos() << std::endl;
      });

    if (argc < 2) return 1;
    // Feed as many file as you want by calling multiple times:
    f2kmers.read_fastx(argv[1]); // The callbacks are called at this point
    /* Support manual sequence feeding too
        f2kmers.next_chrom() // Signal new sequence
        f2kmers.feed("ATGCTGGCTGCTATTC"); // Feed data
        f2kmers.feed("GC"); // Just two new kmers
    */

    std::cout << "Seen " << id_seq << " sequences\n";
}
