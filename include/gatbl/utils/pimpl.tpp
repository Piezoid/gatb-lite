#ifndef GATBL_PIMPL_TPP
#define GATBL_PIMPL_TPP

#include "type_erasure.hpp"

namespace gatbl {

template<typename T>
template<typename... Args>
pimpl<T>::pimpl(Args&&... args)
  : _m(new T{std::forward<Args>(args)...})
{}

template<typename T>
pimpl<T>::pimpl(T* p)
  : _m(p)
{}

template<typename T> pimpl<T>::~pimpl() {}

} // namespace gatbl

#endif // GATBL_PIMPL_TPP
