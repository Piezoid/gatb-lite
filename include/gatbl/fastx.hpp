#ifndef FASTX_HPP
#define FASTX_HPP

#include "gatbl/utils/iterator_pair.hpp"
#include "gatbl/utils/condition_check.hpp"

namespace gatbl {

struct parse_error : std::runtime_error
{
    explicit parse_error(const std::string& what)
      : runtime_error("Parse error: " + what)
    {}
};

template<typename Record> class sequence_iterator : private utils::condition_check
{
    using char_iterator        = typename Record::char_iterator;
    using char_iterator_traits = std::iterator_traits<char_iterator>;

  public:
    using value_type        = Record;
    using reference         = Record&;
    using iterator_category = typename char_iterator_traits::iterator_category;
    using difference_type   = typename char_iterator_traits::difference_type;
    using pointer           = typename char_iterator_traits::pointer;

    /// Initialize an iterator for a sub range of text data: we must first sync to the first header.
    sequence_iterator(const char_iterator& begin, const char_iterator& end)
      : _record(begin)
    {
        _record.sync(end);
    }

    /// Initialize an iterator for a complete file: the iterator must point to a sequence.
    sequence_iterator(char_iterator it)
      : _record(it)
    {
        if (unlikely(*it != value_type::header_char)) { throw parse_error("No header found at the start of the file"); }
    }

    operator char_iterator() { return _record.begin(); }

    sequence_iterator(const sequence_iterator&) = default;
    sequence_iterator(sequence_iterator&&)      = default;
    sequence_iterator& operator=(const sequence_iterator&) = default;
    sequence_iterator& operator=(sequence_iterator&&) = default;

    /// Here is the trick: we parse when checking for completion
    /// This allows to have the sentinel in hand for parsing
    bool operator<(const char_iterator& end) const
    {
        bool parseok = _record.parse(end);
        condition_check::checked();
        return parseok;
    }
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

    bool operator<(const sequence_iterator& other) const { return operator<(other._record.begin()); }

    bool operator!=(const char_iterator& other) const { return operator<(other); }
    bool operator!=(const sequence_iterator& other) const { return operator<(other); }

  private:
    mutable value_type _record;
};

template<typename Record>
using sequence_range = iterator_pair<sequence_iterator<Record>, typename Record::char_iterator>;

template<typename Record>
class default_splitter<sequence_range<Record>> : public default_splitter<iterator_pair<typename Record::char_iterator>>
{
    using base = default_splitter<iterator_pair<typename Record::char_iterator>>;

  public:
    using base::base;
    using typename base::iterator;
    using typename base::sentinel;
    iterator split(iterator beg, sentinel end)
    {
        size_t dist = distance(beg, end);
        if (dist > base::_size) {
            advance(beg, dist / 2);
            Record rec{beg};
            rec.sync(end);
            return begin(rec);
        } else {
            return end;
        }
    }

    sequence_range<Record> make_range(iterator beg, sentinel end) { return {beg, end}; }
};

template<typename CharIt = const char*> struct line_record : iterator_pair<CharIt>
{
    using char_iterator = CharIt;
    using substring_t   = iterator_pair<char_iterator>;

    line_record(CharIt it)
      : substring_t(it, it)
    {
        assume(it, "no input data");
    }

    void check_sync() const {};

    void sync(CharIt end)
    {
        CharIt it = find(substring_t::begin(), end, '\n');
        if (it != end && ++it != end) { // skip '\n'
            substring_t::operator=({it, find(it, end, '\n')});
        } else {
            substring_t::operator=({end, end});
        }
    }

    bool parse(CharIt end)
    {
        CharIt it = substring_t::begin();
        if (it != end && (unlikely(*it != '\n') || ++it != end)) { // skip '\n'
            substring_t::operator=({it, find(it, end, '\n')});
            return substring_t::end() != end;
        } else {
            return false;
        }
    }

    void next() { substring_t::operator=({substring_t::end(), substring_t::end()}); }
};

template<typename CharIt = const char*, typename std::iterator_traits<CharIt>::value_type HeaderChar = '>'>
class fasta_record
{
    using char_iterator_traits = std::iterator_traits<CharIt>;
    static_assert(
      std::is_same<typename char_iterator_traits::iterator_category, std::random_access_iterator_tag>::value,
      "Random acces iterator required");

  public:
    using char_iterator = CharIt;
    using substring_t   = iterator_pair<char_iterator>;
    static_assert(std::is_same<char, typename char_iterator_traits::value_type>::value,
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
        return {_sequence_header_end + 1 /* '\n' */, _sequence_end};
    }

    char_iterator end() const { return _sequence_end; }

    fasta_record(char_iterator it)
      : _data(it)
    {
        assume(_data, "no input data");
    }

    void sync(const char_iterator& end)
    {
        char_iterator p = _data;
        do {
            if (unlikely(p >= end)) {
                _data = end;
                return;
            }
            if (*p == header_char) {
                _data = p;
                return;
            }
            p = find(p, end, '\n') + 1;
        } while (true);
    }

    bool parse(const char_iterator& end)
    {
        if (likely(_parse(end) < end)) { return true; }

        _sequence_header_end = this->begin();
        _sequence_end        = this->begin();
        throw_parse_error();
        return false;
    }

    void next() { this->_data = sequence().end() + 1 /* '\n' */; }

  protected:
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
    char_iterator _parse_header(const char_iterator& end)
    {
        char_iterator p = begin();

        // Sequence header
        if (unlikely(p + 1 >= end)) return end;
        if (unlikely(*(p++) != header_char)) {
            set_error("Expected header");
            return end;
        }
        char_iterator newline = find(p, end, '\n');
        _sequence_header_end  = newline;
        assert(!contains(p, newline, '\n'), "newline found in header");

        return newline;
    }

    /// Returns an iterator to the last parsed char
    char_iterator _parse(const char_iterator& end)
    {
        // Sequence Header
        char_iterator newline = _parse_header(end);
        // Sequence
        char_iterator p = newline + 1; // skip '\n'
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

    void sync(const char_iterator& end)
    {
        base::sync(end);
        char_iterator p = base::begin();
        // Checks if the next line also starts with an '@'
        p = find(p, end, '\n') + 1;
        if (unlikely(p >= end)) {
            base::_data = end;
            return;
        }
        if (*p == base::header_char) {
            // In which case we were on a quality starting with '@' and the header is the next line:
            base::_data = p;
        }
    }

    bool parse(const char_iterator end)
    {
        if (likely(_parse(end) < end)) { return true; }

        this->_sequence_header_end = this->begin();
        this->_sequence_end        = this->begin();
        this->_quality_start       = this->begin();
        base::throw_parse_error();
        return false;
    }

    void next() { this->_data = quality().end() + 1 /* '\n' */; }

  private:
    char_iterator _parse(const char_iterator end)
    {
        char_iterator p = base::_parse(end);
        if (unlikely(p >= end)) return end;
        p++; // '\n'

        // Quality header
        if (unlikely(p + 2 >= end)) return end;
        if (unlikely(*(p++) != '+')) {
            base::set_error("Expected quality header");
            return end;
        }
        char_iterator newline = *(p + 1) == '\n' ? p + 1 : find(p, end, '\n');
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
