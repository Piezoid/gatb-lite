#ifndef GATBL_REVERSE_RANGE_HPP
#define GATBL_REVERSE_RANGE_HPP

#include "gatbl/utils/ranges.hpp"

namespace gatbl {

template<typename Repr>
class reverse_range
  : public view_facade<reverse_range<Repr>,
                       std::reverse_iterator<sentinel_t<Repr>>,
                       std::reverse_iterator<iterator_t<Repr>>>
{
    using base = view_facade<reverse_range<Repr>,
                             std::reverse_iterator<sentinel_t<Repr>>,
                             std::reverse_iterator<iterator_t<Repr>>>;

    Repr _data;

  public:
    using typename base::iterator;
    using typename base::sentinel;
    using reverse_sentinel = sentinel_t<Repr>;
    using reverse_iterator = iterator_t<Repr>;

    template<typename... Args, typename = enable_if_t<std::is_constructible<Repr, Args...>::value>>
    reverse_range(Args&&... args)
      : _data(std::forward<Args>(args)...)
    {}

    iterator begin() const { return iterator{_data.end()}; }
    sentinel end() const { return sentinel{_data.begin()}; }

    reverse_iterator rbegin() const { return _data.begin(); }
    reverse_sentinel rend() const { return _data.end(); }
};

} // namespace gatbl

#endif // GATBL_REVERSE_RANGE_HPP
