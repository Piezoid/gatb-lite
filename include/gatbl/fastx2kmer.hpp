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
    using window_t = Window;
    using typename Window::value_type;
    using Window::size;

    template<typename W>
    sequence2kmers_base(W&& w)
      : window_t(std::forward<W>(w))
      , _unfilled_nucs(window_t::size())
    {}

    ~sequence2kmers_base() { clear(); }

    void clear()
    {
        empty_window();
        _pos = 0;
    }

    /// Push a sequence string
    /// This can be called multiple times for multiple parts of a single sequence
    template<typename R> void feed(R&& r)
    {
        ksize_t       unfilled_nucs = _unfilled_nucs;
        const ksize_t window_size   = window_t::size();

        for (nuc_t n : dna_ascii_range<R>(std::forward<R>(r))) {
            if (unlikely(n == nuc_t::N)) { // Invalid nucletoide
                if (unfilled_nucs == 0) {  // We had a valid kmer
                    derived().on_run_end();
                }
                unfilled_nucs = window_size;

                ++_pos;
                continue;
            } else if (unlikely(unfilled_nucs > 0)) { // Window not full
                Window::unchecked_push_back(n);
                if (--unfilled_nucs > 0) {
                    ++_pos;
                    continue;
                } else { // The window is now full
                    Window::mask(true);
                    Window::check();
                }
            } else { // Slide a full window
                Window::push_back(n);
            }
            ++_pos;
            derived().on_kmer();
        }
        _unfilled_nucs = unfilled_nucs;
    }

    /// Advance position by the length of the sequence
    template<typename R> void skip(R&& r)
    {
        assert(_unfilled_nucs == window_t::size(), "window is not empty, call next_chrom() first");
        using gatbl::size;
        _pos += size(r);
    }

    /// Signal the start of new a chromosome/read, returns true if on_chrom() choose to keep it
    bool start()
    {
        empty_window();
        return derived().on_chrom();
    }

    // Adds a full sequence
    template<typename R> void sequence(R&& r)
    {
        if (start()) {
            feed(std::forward<R>(r));
        } else {
            skip(std::forward<R>(r));
        }
    }

    /// Read a fasta or a fastx and feeds all the reads
    void read_fastx(const std::string& fin)
    {
        gatbl::sys::file_descriptor fd(fin);
        auto                        content = fd.mmap<const char>();
        content.advise_hugepage();

        if (hasEnding(fin, ".fq") || hasEnding(fin, ".fastq")) {
            RANGES_FOR(auto& rec, sequence_range<fastq_record<>>(content)) { sequence(rec.sequence()); }
        } else if (hasEnding(fin, ".fa") || hasEnding(fin, ".fasta")) {
            bool dontskip = false;
            RANGES_FOR(auto& line, sequence_range<line_record<>>(content))
            {
                if (unlikely(gatbl::size(line) == 0)) continue;
                auto it = begin(line);
                if (*it != '>') {
                    if (dontskip) {
                        feed(line);
                    } else {
                        skip(line);
                    }
                } else {
                    ++it;
                    dontskip = start();
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

    const window_t& get_window() const { return static_cast<const window_t&>(*this); }

  protected:
    void empty_window()
    {
        if (_unfilled_nucs == 0) { // We had a valid kmer
            derived().on_run_end();
        }
        _unfilled_nucs = window_t::size();
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
    bool on_chrom() { return _on_chrom(as_ref()); }

    const sequence2kmers& as_ref() const { return *this; }

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

    // We don't want functors to make accidental copies (infinite recursion)
    sequence2kmers(const sequence2kmers&) = delete;
    sequence2kmers& operator=(const sequence2kmers&) = delete;
    sequence2kmers(sequence2kmers&&)                 = default;
    sequence2kmers& operator=(sequence2kmers&&) = default;
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
