#ifndef GATBL_FASTX2KMER_HPP
#define GATBL_FASTX2KMER_HPP

#include <gatbl/utils/nop_functor.hpp>
#include <gatbl/utils/string.hpp>

#include <gatbl/sys/file.hpp>
#include <gatbl/kmer.hpp>
#include <gatbl/fastx.hpp>

namespace gatbl {
/// Convert sequence to kmers, with capabilities to avoid invalid chars and multiline FASTA
template<typename KmerT = kmer_t, typename OnKmer = nop_functor, typename OnRunEnd = nop_functor> struct sequence2kmers
{
    using kmer_t = KmerT;
    using arg_t  = positioned<sized_kmer<kmer_t>>;

    sequence2kmers(ksize_t k, OnKmer&& on_kmer = {}, OnRunEnd&& on_run_end = {})

      : _window(k)
      , _window_unfilled_nucs(_window.size())
      , _on_kmer(std::move(on_kmer))
      , _on_run_end(std::move(on_run_end))
    {}

    // Push a sequence string

    template<typename R> void feed(const R& r)
    {
        dna_ascii_range<R> seq(r);
        // debug_op(std::cerr << "s: " << r << endl);

        auto       it   = begin(seq);
        const auto last = end(seq);
        if (unlikely(_window_unfilled_nucs != 0)) it = fill(it, last);

        for (; it != last;) {
            while (unlikely(*it == nuc_t::N)) {
                empty_window();
                it = fill(it, last);
                if (it == last) return;
            }

            _window.push_back(*it);
            ++_pos;
            ++it;

            _on_kmer(arg_t{_window.forward(), _pos});
        }
    }

    // Get the position of the next nucleotide
    size_t get_next_pos() const { return _pos; }

    // Signal the start of a chromosome/FASTA/Q sequence
    void next_chrom()
    {
        _chrom_starts.push_back(_pos);
        empty_window();
    }

    void read_fastx(const std::string fin)
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
                if (unlikely(size(line) == 0)) continue;
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

    // std::vector<size_t>&&      get_chrom_starts() && { return std::move(_chrom_starts); }
    std::vector<size_t>&       get_chrom_starts() { return _chrom_starts; }
    const std::vector<size_t>& get_chrom_starts() const { return _chrom_starts; }

  protected:
    void empty_window()
    {
        if (_window_unfilled_nucs == 0) { // We had a kmer
            _on_run_end(arg_t{_window.forward(), _pos});
        }
        _window_unfilled_nucs = _window.size();
    }

    /// Fill the kmer window at the begining of a line or after Ns.
    /// Filling can be done in mulitple call (eg when Ns span multiple lines)
    template<typename It, typename S> It fill(It it, const S last)
    {
        for (; _window_unfilled_nucs > 0 && it != last; ++it, ++_pos) {
            if (*it == nuc_t::N) { // Ns likely to follow N
                empty_window();
            } else {
                _window.unchecked_push_back(*it);
                --_window_unfilled_nucs;
            }
        }
        if (_window_unfilled_nucs == 0) {
            _window.mask(true);
            _window.check();
            _on_kmer(arg_t{_window.forward(), _pos});
        }
        return it;
    }

    std::vector<size_t> _chrom_starts;
    kmer_window<kmer_t> _window;
    size_t              _pos = 0;
    ksize_t             _window_unfilled_nucs;
    OnKmer              _on_kmer;
    OnRunEnd            _on_run_end;
};

template<typename KmerT = kmer_t, typename OnKmer = nop_functor, typename OnRunEnd = nop_functor>
sequence2kmers<KmerT, remove_reference_t<OnKmer>, remove_reference_t<OnRunEnd>>
make_sequence2kmers(ksize_t k, OnKmer&& on_kmer = {}, OnRunEnd&& on_run_end = {})
{
    return {k, std::forward<OnKmer>(on_kmer), std::forward<OnRunEnd>(on_run_end)};
}

}

#endif // GATBL_FASTX2KMER_HPP
