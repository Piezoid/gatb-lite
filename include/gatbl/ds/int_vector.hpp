#ifndef GATBL_INT_VECTOR_HPP
#define GATBL_INT_VECTOR_HPP

#include "gatbl/utils/ranges.hpp"
#include "gatbl/sys/bits.hpp"
#include "gatbl/utils/dna.hpp" // For underlying type and dna enum

namespace gatbl {

namespace details {

using bits::bitsize_t;

/// Represents the width in bits of integers, either dynamically or by a type constant
template<typename WidthT = bitsize_t> class bitwidth
{
    static constexpr size_t width_bits = CHAR_BIT * sizeof(WidthT);
    static constexpr size_t width_mask = (size_t(1) << width_bits) - 1;

    WidthT _width;

  public:
    constexpr bitwidth(WidthT width = {}) noexcept
      : _width(width)
    {}
    template<typename T> static constexpr WidthT forvalue(T v) noexcept
    {
        return {WidthT(bits::ilog2p1(to_underlying(v)))};
    }
    WidthT constexpr width() const noexcept { return _width; }

    // For serialization
    constexpr size_t       pack_size_and_width(size_t size) const { return (size << width_bits) | _width; }
    CPP14_CONSTEXPR size_t unpack_size_and_width(size_t size)
    {
        _width = WidthT(size & width_mask);
        return size >> width_bits;
    }
};

/// Statically represents the width in bits of integers. Will occupy 0 bytes when used as a base class.
template<typename WidthT, WidthT N> class bitwidth<std::integral_constant<WidthT, N>>
{
  public:
    constexpr bitwidth(std::integral_constant<WidthT, N> = {}) noexcept {}
    template<typename T> static CPP14_CONSTEXPR std::integral_constant<WidthT, N> forvalue(T v) noexcept
    {
        assume(to_underlying(v) < (underlying_type<T>(1) << N), "value out of bounds");
        return {};
    }
    constexpr std::integral_constant<WidthT, N> width() const noexcept { return {}; }

    // For serialization
    constexpr size_t pack_size_and_width(size_t size) const { return size; }
    constexpr size_t unpack_size_and_width(size_t size) const { return size; }
};

/// Specifies the representation for integers packed in memory
/// Representatively isomorphic to bitwidth<WidthT>, but with helper methods
/// FIXME: this version relies heavilly on little-endian and unaligned word access: the implementation for other
/// architectures is present but untested and probably ineficient.
template<typename T      = size_t,              /// Value type.
         typename WidthT = bitsize_t,           /// Type representing the integer width in bits,
                                                /// can be std::integral_constant<uint8_t, N>.
         typename WordT = size_t>               /// Word type loaded from memory (char-aligned!).
struct packedint_spec : public bitwidth<WidthT> // bitwidth specification, EBOptimized
{
    using base    = bitwidth<WidthT>;
    using width_t = WidthT;
    using word_t  = WordT;
    using type    = T;

    static constexpr bitsize_t word_bits  = sizeof(word_t) * CHAR_BIT;
    static constexpr bitsize_t offset_max = CHAR_BIT - 1;

    CPP14_CONSTEXPR packedint_spec(WidthT width = {}) noexcept
      : base(width)
    {
        check_offset(offset_max);
    }

    CPP14_CONSTEXPR void check_offset(bitsize_t offset) const noexcept
    {
        assume(offset <= offset_max && base::width() + offset <= word_bits, "word overflow");
    }

    constexpr word_t mask() const noexcept { return word_t(word_t(-1) >> (word_bits - base::width())); }

    CPP14_CONSTEXPR type load(const word_t* ptr, bitsize_t offset) const noexcept
    {
        check_offset(offset);

        auto word = bits::load_le(ptr);
        return type((word >> offset) & mask());
        //        word <<= bitsize_t(word_bits - width() - offset);
        //        word >>= bitsize_t(word_bits - width());
        //        return type(word);
    }

    CPP14_CONSTEXPR void store(word_t* ptr, bitsize_t offset, type _value) const noexcept
    {
        check_offset(offset);

        word_t       word  = bits::load_le(ptr);
        const word_t m     = mask();
        const word_t value = static_cast<word_t>(_value);
        assume(value <= m, "value out of bounds");
        word = word_t(word & ~(m << offset));
        word = word_t(word | (value << offset));
        bits::store_le(ptr, word);
    }

    /// Returns the number of bits required to store `n` integers of width `width`
    constexpr size_t words_required(size_t size) const noexcept
    {
        // const size_t bits_used       = (size - 1) * base::width();
        // const size_t max_byte_offset = bits_used / CHAR_BIT;
        // const size_t required_bytes  = max_byte_offset + sizeof(word_t);
        // const size_t required_words  = (required_bytes + (sizeof(word_t) - 1)) / sizeof(word_t);
        return size > 0 ? ((size - 1) * base::width() + (2 * word_bits - CHAR_BIT)) / word_bits : 0;
    }
};

using dynamic_packed_int_spec               = packedint_spec<>;
template<bitsize_t N> using static_bitwidth = std::integral_constant<bitsize_t, N>;
using packed_twobit_spec                    = packedint_spec<nuc_t, static_bitwidth<2>>;
using bit_spec                              = packedint_spec<bool, static_bitwidth<1>>;

/// Base namespace declaring types for packed int containers
template<typename Spec = dynamic_packed_int_spec> struct packedint_types
{
    struct iterator;
    struct reference;
    struct const_iterator;
    using const_reference = typename Spec::type;

  protected:
    using bitptr_t = std::ptrdiff_t; // An absolute bit address. A signed type allow arithmetic right shift to restore a
                                     // cannonical x86_64 pointer

    /// Common CRTP facade
    template<typename D, typename Ref>
    class iterator_base
      : public iterator_facade<D, std::random_access_iterator_tag, typename Spec::type, Ref, void, bitptr_t>
      , protected Spec // bitwidth specification, EBOptimized
    {
      protected:
        using base = iterator_facade<D, std::random_access_iterator_tag, typename Spec::type, Ref, void, bitptr_t>;

        bitptr_t addr = 0;


      public:
        constexpr iterator_base() noexcept = default;

        constexpr iterator_base(bitptr_t bitptr, Spec spec) noexcept
          : Spec(spec)
          , addr(bitptr)
        {}

        CPP14_CONSTEXPR iterator_base(const typename Spec::word_t* p, Spec w, size_t i = 0) noexcept
          : iterator_base((reinterpret_cast<bitptr_t>(p) << 3) + bitptr_t(i * w.width()), w)
        {
            assume(Spec::width() > 0, "Can't iterate integers of 0 width");
        }

        CPP14_CONSTEXPR D& operator+=(bitptr_t x) noexcept
        {
            this->addr += x * Spec::width();
            return static_cast<D&>(*this);
        }

        friend CPP14_CONSTEXPR bitptr_t compare(const D& lhs, const D& rhs) noexcept
        {
            bitptr_t d = lhs.addr - rhs.addr;
            assume(lhs.width() == rhs.width() && d % lhs.width() == 0, "unaligned addresses compared");
            return d;
        }

        using base::             operator-;
        CPP14_CONSTEXPR bitptr_t operator-(const D& other) const noexcept
        {
            return compare(static_cast<const D&>(*this), other) / Spec::width();
        }
    };

    using word_t = typename Spec::word_t;

  public:
    struct const_iterator : iterator_base<const_iterator, const_reference>
    {
        /// Conversion from non-const iterator
        constexpr const_iterator(const iterator& it) noexcept
          : iterator_base<const_iterator, const_reference>(it.addr, it)
        {}

        using iterator_base<const_iterator, const_reference>::iterator_base;
        constexpr const_reference operator*() const noexcept
        {
            return Spec::load(reinterpret_cast<const word_t*>(this->addr >> 3), bitsize_t(this->addr & 7u));
        }
    };

    struct iterator : iterator_base<iterator, reference>
    {
        friend struct const_iterator;
        using iterator_base<iterator, reference>::iterator_base;
        constexpr reference operator*() const noexcept
        {
            return {reinterpret_cast<word_t*>(this->addr >> 3), bitsize_t(this->addr & 7u), *this};
        }
    };

    struct reference : private Spec // bitwidth specification, EBOptimized
    {
        word_t*   _ptr;
        bitsize_t _offset;

      public:
        CPP14_CONSTEXPR reference(word_t* ptr, bitsize_t offset, Spec width = {}) noexcept
          : Spec(width)
          , _ptr(ptr)
          , _offset(offset)
        {
            Spec::check_offset(offset);
        }
        using type = typename Spec::type;

        CPP14_CONSTEXPR operator type() const noexcept { return Spec::load(_ptr, _offset); }

        CPP14_CONSTEXPR const reference& operator=(type value) const noexcept
        {
            Spec::store(_ptr, _offset, value);
            return *this;
        }

        CPP14_CONSTEXPR const reference& operator=(const reference& from) const noexcept
        {
            Spec::store(_ptr, _offset, Spec::load(from._ptr, from._offset));
            return *this;
        }

        friend CPP14_CONSTEXPR void swap(reference ref_x, reference ref_y) noexcept
        { // FIXME: masked XOR swap may be faster
            type x = ref_x.load(ref_x._ptr, ref_x._offset);
            type y = ref_y.load(ref_y._ptr, ref_y._offset);
            ref_x.store(ref_x._ptr, ref_x._offset, y);
            ref_y.store(ref_y._ptr, ref_y._offset, x);
        }
    };
};

}

template<typename Spec = details::dynamic_packed_int_spec>
class int_vector
  : public view_facade<int_vector<Spec>,
                       typename details::packedint_types<Spec>::iterator,
                       typename details::packedint_types<Spec>::const_iterator>
  , protected Spec // bitwidth specification, EBOptimized
{
  protected:
    using word_t     = typename Spec::word_t;
    using types_base = details::packedint_types<Spec>;

    std::unique_ptr<word_t[]> _data{};
    size_t                    _nwords = 0;
    size_t                    _size   = 0;

  public:
    using value_type = typename Spec::type;
    /// The type used to represent the integer bitwidth, can be implicitly converted from integer, or default
    /// constructed for static width
    using typename Spec::width_t;

    constexpr int_vector() noexcept                   = default;
    CPP14_CONSTEXPR int_vector(int_vector&&) noexcept = default;
    CPP14_CONSTEXPR int_vector& operator=(int_vector&&) noexcept = default;

    /// Allocate a vector of size integers, no initialization if performed
    CPP14_CONSTEXPR int_vector(size_t size, width_t width = {})
      : Spec(width)
    {
        resize(size, false, false);
    }

    /// Initialize with a fill value
    CPP14_CONSTEXPR int_vector(size_t size, width_t width, value_type fill_value)
      : int_vector(size, width)
    {
        if (size_t(fill_value) == 0u) {
            std::fill(_data.get(), _data.get() + _nwords, 0);
        } else {
            std::fill(this->begin(), this->end(), fill_value);
        }
    }

    CPP14_CONSTEXPR int_vector(const int_vector& from)
      : Spec(from)
    {
        resize(from._size, false);
        const auto* from_ptr = from._data.get();
        std::copy(from_ptr, from_ptr + _nwords, _data.get());
    }

    CPP14_CONSTEXPR int_vector& operator=(const int_vector& from)
    {
        static_cast<Spec*>(this)->operator=(from);
        resize(from._size, false);
        const auto* from_ptr = from._data.get();
        std::copy(from_ptr, from_ptr + _nwords, _data.get());
        return *this;
    }

    /// Initialize a vector from a range
    /// The maximum value will set the integer bitwidth or be checked against the static bitwidth
    template<typename R, typename = decltype(concepts::RangeOf<const typename Spec::type, R>)>
    CPP14_CONSTEXPR int_vector(const R& r, typename R::value_type max_value)
      : int_vector(r.size(), Spec::forvalue(max_value))
    {
        using std::begin;
        using std::end;
        std::copy(begin(r), end(r), this->begin());
    }

    /// Initialize a vector from a range, maximal value is obtained from std::max_element
    template<typename R>
    CPP14_CONSTEXPR explicit int_vector(const R& r)
      : int_vector(r, *std::max_element(r.begin(), r.end()))
    {}

    constexpr size_t size() const noexcept { return _size; }
    constexpr bool   empty() const noexcept { return _size == 0; }
    using Spec::width;

    CPP14_CONSTEXPR void resize(size_t size, bool keep_data = true, bool shrink = false)
    {
        _size           = size;
        auto new_nwords = Spec::words_required(size);
        if ((new_nwords > _nwords) || (shrink && new_nwords < _nwords)) { // realloc
            auto   new_data = std::unique_ptr<word_t[]>(new word_t[new_nwords]);
            size_t copysize = new_nwords > _nwords ? _nwords : new_nwords;
            if (keep_data) std::copy(_data.get(), _data.get() + copysize, new_data.get());
            _nwords = new_nwords;
            _data   = std::move(new_data);
        }
    }

    CPP14_CONSTEXPR void shrink_to_fit(bool keep_data = true) { resize(_size, keep_data, true); }

    constexpr size_t size_in_bytes() const { return _nwords * sizeof(word_t); }

    CPP14_CONSTEXPR typename types_base::const_iterator begin() const noexcept { return {_data.get(), *this}; }
    CPP14_CONSTEXPR typename types_base::const_iterator end() const noexcept { return {_data.get(), *this, _size}; }
    CPP14_CONSTEXPR typename types_base::iterator       begin() noexcept { return {_data.get(), *this}; }
    CPP14_CONSTEXPR typename types_base::iterator       end() noexcept { return {_data.get(), *this, _size}; }

#ifdef GATBL_RANGE_IO_HPP
    template<typename O> friend O& write(O& out, const int_vector& v)
    {
        size_t header = v.Spec::pack_size_and_width(v._size);
        write(out, header);
        write(out, span<const word_t>(v._data.get(), v.words_required(v._size)));
        return out;
    }

    template<typename I> friend I& read(I& in, int_vector& v)
    {
        size_t header;
        read(in, header);
        v.resize(v.Spec::unpack_size_and_width(header), false, false);
        read(in, span<word_t>(v._data.get(), v.words_required(v._size)));
        return in;
    }
#endif
};

using twobit_dna_vector = int_vector<details::packed_twobit_spec>;
using bit_vector        = int_vector<details::bit_spec>;

} // namespace gatbl;

#endif // GATBL_INT_VECTOR_HPP
