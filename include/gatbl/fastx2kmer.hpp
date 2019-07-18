#ifndef GATBL_FASTX2KMER_HPP
#define GATBL_FASTX2KMER_HPP

#include <gatbl/utils/nop_functor.hpp>
#include <gatbl/utils/string.hpp>

#include <gatbl/sys/file.hpp>
#include <gatbl/kmer.hpp>
#include <gatbl/fastx.hpp>

namespace gatbl {

namespace details {

/// A high level fastx to kmers
/// It handles fastq, genomic multiline fasta parsing and the presence of invalid nucleotide in the range
template<typename Derived,                     // CRTP Derived type
         typename Window = kmer_window<kmer_t> // Window used to accumulate nucleotides
         >
struct sequence2kmers_base : private Window
{
    template<typename W>
    sequence2kmers_base(W&& w)
      : Window(std::forward<W>(w))
      , _unfilled_nucs(Window::size())
    {}

    using window_t = Window;
    using typename Window::value_type;
    using Window::size;

    /// Push a sequence string
    /// This can be called multiple times for multiple parts of a single sequence
    template<typename R> void feed(R&& r)
    {
        dna_ascii_range<R> seq(r);

        auto       it   = begin(seq);
        const auto last = end(seq);
        if (unlikely(_unfilled_nucs != 0)) it = fill(it, last);

        for (; it != last;) {
            while (unlikely(*it == nuc_t::N)) {
                empty_window();
                it = fill(it, last);
                if (it == last) return;
            }

            Window::push_back(*it);
            ++it;
            ++_pos;

            derived().on_kmer();
        }
    }

    /// Signal the start of new a chromosome/read
    void next_chrom()
    {
        derived().on_chrom();
        empty_window();
    }

    /// Read a fasta or a fastx and feeds all the reads
    void read_fastx(const std::string& fin)
    {
        gatbl::sys::file_descriptor fd(fin);
        auto                        content = fd.mmap<const char>();
        content.advise_hugepage();

        if (hasEnding(fin, ".fq") || hasEnding(fin, ".fastq")) {
            RANGES_FOR(auto& rec, sequence_range<fastq_record<>>(content))
            {
                next_chrom();
                feed(rec.sequence());
            }
        } else if (hasEnding(fin, ".fa") || hasEnding(fin, ".fasta")) {
            RANGES_FOR(auto& line, sequence_range<line_record<>>(content))
            {
                if (unlikely(gatbl::size(line) == 0)) continue;
                auto it = begin(line);
                if (*it != '>') {
                    feed(line);
                } else {
                    ++it;
                    next_chrom();
                }
            }
        } else {
            throw std::runtime_error("unsupported file format");
        }
    }

    /// Get the position of the next nucleotide to be added to the window, in global coordinates (for all
    /// reads/chrom)
    size_t get_next_pos() const { return _pos; }

    /// Get the position at the left of the first kmer (enough sequence must be fed before)
    size_t get_pos() const
    {
        assert(_pos >= Window::size(), "Invalid position, is the window filled ?");
        return _pos - Window::size();
    }

    const Window& get_window() const { return static_cast<const Window&>(*this); }

  protected:
    void empty_window()
    {
        if (_unfilled_nucs == 0) { // We had a valid kmer
            derived().on_run_end();
        }
        _unfilled_nucs = Window::size();
    }

    /// Fill the kmer window at the begining of a line or after Ns.
    /// Filling can be done in multople call (eg. when Ns span multiple lines)
    template<typename It, typename S> It fill(It it, const S last)
    {
        for (; _unfilled_nucs > 0 && it != last; ++it, ++_pos) {
            if (*it == nuc_t::N) { // Ns likely to follow N
                empty_window();
            } else {
                Window::unchecked_push_back(*it);
                --_unfilled_nucs;
            }
        }
        if (_unfilled_nucs == 0) {
            Window::mask(true);
            Window::check();
            derived().on_kmer();
        }
        return it;
    }

  private:
    Derived& derived() { return static_cast<Derived&>(*this); }
    Derived& derived() const { return static_cast<Derived&>(*this); }

    size_t  _pos = 0;
    ksize_t _unfilled_nucs;
};

} // namespace details

/// Convert sequence to kmers, with capabilities to avoid invalid chars and multiline FASTA
template<typename Window   = kmer_window<kmer_t>,
         typename OnKmer   = nop_functor, // Lambda called with windows when a new kmer is seen
         typename OnChrom  = nop_functor, // Lambda called with the last position before another sequence
         typename OnRunEnd = nop_functor  // Lambda called with the last kmer just before
                                          // the end of a sequence or a invalide nucleotide
         >
class sequence2kmers : public details::sequence2kmers_base<sequence2kmers<Window, OnKmer, OnChrom, OnRunEnd>, Window>
{
    using base = details::sequence2kmers_base<sequence2kmers<Window, OnKmer, OnChrom, OnRunEnd>, Window>;

    friend base;
    void on_kmer() { _on_kmer(as_ref()); }
    void on_run_end() { _on_run_end(as_ref()); }
    void on_chrom() { _on_chrom(as_ref()); }

    OnKmer   _on_kmer;
    OnChrom  _on_chrom;
    OnRunEnd _on_run_end;

  public:
    template<typename W>
    sequence2kmers(W&& w, OnKmer&& on_kmer = {}, OnChrom&& on_chrom = {}, OnRunEnd&& on_run_end = {})
      : base(w)
      , _on_kmer(std::forward<OnKmer>(on_kmer))
      , _on_chrom(std::forward<OnChrom>(on_chrom))
      , _on_run_end(std::forward<OnRunEnd>(on_run_end))
    {}

    const sequence2kmers& as_ref() const { return *this; }
};

template<typename Window   = kmer_window<kmer_t>,
         typename OnKmer   = nop_functor,
         typename OnChrom  = nop_functor,
         typename OnRunEnd = nop_functor>
sequence2kmers<Window, OnKmer, OnChrom, OnRunEnd>
make_sequence2kmers(Window&& w, OnKmer&& on_kmer = {}, OnChrom&& on_chrom = {}, OnRunEnd&& on_run_end = {})
{
    return {std::forward<Window>(w),
            std::forward<OnKmer>(on_kmer),
            std::forward<OnChrom>(on_chrom),
            std::forward<OnRunEnd>(on_run_end)};
}
}

#endif // GATBL_FASTX2KMER_HPP
