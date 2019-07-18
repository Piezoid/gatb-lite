
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

      // Callback for new kmer
      [&](auto& f2kmers) {
          sized_kmer<kmer_t> kmer = f2kmers.get_window().forward();
          std::cout << kmer // sized_kmer is printable since it know it real size
                    << " " << f2kmers.get_pos() << std::endl;
      },

      // New sequence in fast[aq]
      [&](auto& f2kmers) {
          bool keep = (id_seq % 2) == 1; // only make kmers on odd sequence
          std::cout << "sequence " << id_seq++ << " start:" << f2kmers.get_next_pos() << " keep=" << keep << std::endl;
          return true; // return true to keep the sequence
      },
      // Called before the windows is emptied, when sequence is interrupted by a N or at the last nucleotide of a
      // sequence (after a call to start())
      [&](auto& f2kmers) { std::cout << "end run " << f2kmers.get_pos() << std::endl; });

    if (argc < 2) return 1;

    // Feed as many file as you want by calling multiple times:
    f2kmers.read_fastx(argv[1]); // The callbacks are called at this point

    // Support manual sequence feeding too
    f2kmers.start();                               // Signal new sequence
    f2kmers.feed(std::string("ATGCTGGCTGCTATTC")); // Feed data
    f2kmers.feed(std::string("GC"));               // adds two new kmers to the previous sequence

    std::cout << "Seen " << id_seq << " sequences\n";
}
