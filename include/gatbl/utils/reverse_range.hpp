#ifndef GATBL_REVERSE_RANGE_HPP
#define GATBL_REVERSE_RANGE_HPP

#include "gatbl/utils/ranges.hpp"

namespace gatbl {

template<typename Repr>
class reverse_range
  : public view_facade<reverse_range<Repr>,
                       std::reverse_iterator<sentinel_t<Repr>>,
                       std::reverse_iterator<sentinel_t<const Repr>>,
                       std::reverse_iterator<iterator_t<Repr>>,
                       std::reverse_iterator<iterator_t<const Repr>>>
{
    using base = view_facade<reverse_range<Repr>,
                             std::reverse_iterator<sentinel_t<Repr>>,
                             std::reverse_iterator<sentinel_t<const Repr>>,
                             std::reverse_iterator<iterator_t<Repr>>,
                             std::reverse_iterator<iterator_t<const Repr>>>;

    Repr _data;

  public:
    using reverse_iterator       = iterator_t<Repr>;
    using reverse_sentinel       = sentinel_t<Repr>;
    using reverse_const_iterator = iterator_t<const Repr>;
    using reverse_const_sentinel = sentinel_t<const Repr>;

    template<typename... Args, typename = enable_if_t<std::is_constructible<Repr, Args...>::value>>
    reverse_range(Args&&... args)
      : _data(std::forward<Args>(args)...)
    {}

    typename base::const_iterator begin() const { return {_data.end()}; }
    typename base::const_sentinel end() const { return {_data.begin()}; }

    typename base::iterator begin() { return {_data.end()}; }
    typename base::sentinel end() { return {_data.begin()}; }

    reverse_const_iterator rbegin() const { return _data.begin(); }
    reverse_const_sentinel rend() const { return _data.end(); }

    reverse_iterator rbegin() { return _data.begin(); }
    reverse_sentinel rend() { return _data.end(); }
};

} // namespace gatbl

#endif // GATBL_REVERSE_RANGE_HPP
