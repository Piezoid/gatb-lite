#ifndef FASTX_HPP
#define FASTX_HPP

#include "gatbl/utils/interator_pair.hpp"
#include "gatbl/utils/condition_check.hpp"

namespace gatbl {

struct parse_error : std::runtime_error
{
    explicit parse_error(const std::string& what)
      : runtime_error("Parse error: " + what)
    {}
};

template<typename Record> struct sequence_iterator : private utils::condition_check
{
    using char_iterator = typename Record::char_iterator;
    using range         = iterator_pair<char_iterator>;
    using value_type    = Record;

    /// Initialize an iterator for a sub range of text data: we must first sync to the first header.
    sequence_iterator(const char_iterator& begin, const char_iterator& end)
    {
        char_iterator p = find(begin, end, '\n');

        while (++p < end) {
            if (*p == '@') {
                _record = {p};
                return;
            }
            p = find(p, end, '\n');
        }

        _record = {end};
    }

    sequence_iterator(const sequence_iterator&) = default;
    sequence_iterator(sequence_iterator&&)      = default;
    sequence_iterator& operator=(const sequence_iterator&) = default;
    sequence_iterator& operator=(sequence_iterator&&) = default;

    // This pair of implicit conversions allows to have a seemingly different end iterator type
    sequence_iterator(char_iterator it)
      : _record(it)
    {}
    operator char_iterator() const { return _record.begin(); }

    bool operator<(char_iterator& end) const
    {
        bool parseok = _record.parse(end);
        condition_check::set(parseok);
        return parseok;
    }

    bool operator!=(char_iterator& end) const { return operator<(end); }

    sequence_iterator& operator++()
    {
        condition_check::check();
        condition_check::unchecked();
        _record.next();
        return *this;
    }

    value_type& operator*()
    {
        condition_check::check();
        return _record;
    }

  private:
    mutable value_type _record;
};

} // namespace gatbl

template<typename Record> struct std::iterator_traits<gatbl::sequence_iterator<Record>>
{
    using value_type        = gatbl::sequence_iterator<Record>;
    using pointer           = value_type*;
    using reference         = value_type&;
    using iterator_category = std::forward_iterator_tag;
};

namespace gatbl {

template<typename Record>
inline iterator_pair<sequence_iterator<Record>, typename Record::char_iterator>
seq_record_range(const typename Record::char_iterator& begin, const typename Record::char_iterator& end)
{
    return {begin, end};
}

template<typename Record, typename Range>
inline auto
seq_record_range(const Range& char_range) -> decltype(seq_record_range<Record>(char_range.begin(), char_range.end()))
{
    return seq_record_range<Record>(char_range.begin(), char_range.end());
}

template<typename Record>
inline iterator_pair<sequence_iterator<Record>, typename Record::char_iterator>
seq_record_subrange(const typename Record::char_iterator& begin, const typename Record::char_iterator& end)
{
    return {sequence_iterator<Record>(begin, end), end};
}

template<typename Record, typename Range>
inline auto
seq_record_subrange(const Range& char_range)
  -> decltype(seq_record_subrange<Record>(char_range.begin(), char_range.end()))
{
    return seq_record_subrange<Record>(char_range.begin(), char_range.end());
}

template<typename CharIt = const char*, typename std::iterator_traits<CharIt>::value_type HeaderChar = '>'>
class fasta_record
{
    using char_iterator_trais = std::iterator_traits<CharIt>;
    static_assert(std::is_same<typename char_iterator_trais::iterator_category, std::random_access_iterator_tag>::value,
                  "Random acces iterator required");

  public:
    using char_iterator = CharIt;
    using substring_t   = iterator_pair<char_iterator>;
    static_assert(std::is_same<char, typename char_iterator_trais::value_type>::value,
                  "iterator must yield char values");
    static constexpr char header_char = HeaderChar;

    fasta_record()                    = default;
    fasta_record(const fasta_record&) = default;
    fasta_record(fasta_record&&)      = default;
    fasta_record& operator=(const fasta_record&) = default;
    fasta_record& operator=(fasta_record&&) = default;

    operator bool() const { return _sequence_end != _data; }

    char_iterator begin() const { return _data; }

    const substring_t header() const
    {
        assume(*this, "empty fasta record");
        return {_data + 1 /* HeaderChar */, _sequence_header_end};
    }

    const substring_t sequence() const
    {
        assume(*this, "empty fasta record");
        return {_sequence_header_end + 1 /* "\n" */, _sequence_end};
    }

    char_iterator end() const { return _sequence_end; }

  protected:
    friend sequence_iterator<fasta_record>;

    fasta_record(const char_iterator& it)
      : _data(it)
    {}

    bool parse(const char_iterator& start, const char_iterator& end)
    {
        _data = start;
        if (likely(_parse(end) < end)) { return true; }

        _sequence_header_end = start;
        _sequence_end        = start;
        throw_parse_error();
        return false;
    }

    void set_error(const char* e) noexcept
    {
        assume(error == nullptr, "no error set");
        error = e;
    }

    noinline_fun void throw_parse_error()
    {
        if (unlikely(error != nullptr)) {
            throw parse_error(error);
            error = nullptr;
        }
    }

    /// Returns an iterator to the last parsed char
    char_iterator _parse(const char_iterator& end)
    {
        assume(_data && _data <= end, "no input data");
        const char* p = begin();

        // Sequence header
        if (unlikely(p + 1 >= end)) return end;
        if (unlikely(*(p++) != header_char)) {
            set_error("Expected header");
            return end;
        }
        const char* newline  = find(p, end, '\n');
        _sequence_header_end = newline;
        assert(!contains(p, newline, '\n'), "newline found in header");

        // Sequence
        p = newline + 1; // skip "\n"
        if (unlikely(p >= end)) return end;
        if (unlikely(*p == '\n')) {
            set_error("Empty sequence");
            return end;
        }
        newline       = find(p, end, '\n');
        _sequence_end = newline;
        assert(!contains(p, newline, '\n'), "newline found in sequence");

        return newline;
    }

    const char*   error                = nullptr;
    char_iterator _data                = {};
    char_iterator _sequence_header_end = {};
    char_iterator _sequence_end        = {};
};

template<typename CharIt = const char*> class fastq_record : public fasta_record<CharIt, '@'>
{
    using base = fasta_record<CharIt, '@'>;

  public:
    using base::base;
    using typename base::char_iterator;
    using typename base::substring_t;

    const substring_t quality() const { return {_quality_start, _quality_start + this->sequence().size()}; }

    char_iterator end() const { return quality().end(); }

  protected:
    friend sequence_iterator<fastq_record>;

    bool parse(const char_iterator end)
    {
        if (likely(_parse(end) < end)) { return true; }

        this->_sequence_header_end = this->begin();
        this->_sequence_end        = this->begin();
        this->_quality_start       = this->begin();
        base::throw_parse_error();
        return false;
    }

    void next() { this->_data = quality().end() + 1 /* "\n" */; }

  private:
    char_iterator _parse(const char_iterator end)
    {
        char_iterator p = base::_parse(end);
        if (unlikely(p >= end)) return end;
        p++; // "\n"

        // Quality header
        if (unlikely(p + 1 >= end)) return end;
        if (unlikely(*(p++) != '+')) {
            base::set_error("Expected quality header");
            return end;
        }
        char_iterator newline = find(p, end, '\n'); //, initialized ? p + _quality_header_length : p + 1);
        assert(!contains(p, newline, '\n'), "newline found in quality header");

        // Quality follow...
        p              = newline + 1;
        _quality_start = p;

        // Checks for the newline after quality
        newline = p + this->sequence().size();
        if (unlikely(newline >= end)) return end;
        if (unlikely(*newline != '\n')) {
            base::set_error("Expected newline after quality");
            return end;
        }
        assert(!contains(p, newline, '\n'), "newline found in quality");

        return newline;
    }

    char_iterator _quality_start = {};
};

}

#endif // FASTX_HPP
