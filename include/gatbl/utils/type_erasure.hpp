#ifndef GATBL_TYPE_ERASURE
#define GATBL_TYPE_ERASURE

#include <memory>
#include <type_traits>

#include "concepts.hpp"
#include "../common.hpp"

namespace gatbl {

/// Pimpl idiom wrapper, taken from Herb Sutter.
/// This class can be instantiated while T is an incomplete type.
/// You need to include gatbl/utils/pimpl.tpp in your .cpp and instantiate in the footer:
/// template class gatbl::pimpl<MyClassImpl>;
template<typename T> class pimpl
{
    std::unique_ptr<T> _m;
    T&                 _get() const
    {
        T* p = _m.get();
        assume(p != nullptr, "pimpl: null pointer");
        return *p;
    }

  public:
    template<typename... Args> pimpl(Args&&...);
    pimpl(T*);
    ~pimpl();

    pimpl(pimpl&&) = default;
    pimpl& operator=(pimpl&&) = default;
    // FIXME: support copy, maybe ?

    T*       operator->() { return &_get(); }
    const T* operator->() const { return &_get(); }
    T&       operator*() { return _get(); }
    const T& operator*() const { return _get(); }
};

/// Type erased closure reference
/// Allows to define ABI that accept anonymous closures/functor callbacks
/// Beware: no lifetime nor RAII management is done, the closure must remain allocated untill the invocation is made
template<typename Sign> class erased_closure_ref;

// FIXME: constness progpagation: variant for reference to const closures
template<typename Ret, typename... Args> class erased_closure_ref<Ret(Args...)>
{
    Ret (*_fun_ptr)(void*, Args...);
    void* _closure_ptr;

  public:
    erased_closure_ref()                          = delete;
    erased_closure_ref(const erased_closure_ref&) = default;
    erased_closure_ref(erased_closure_ref&)       = default;
    erased_closure_ref& operator=(const erased_closure_ref&) = default;
    erased_closure_ref& operator=(erased_closure_ref&&) = default;

    template<typename F, typename = decltype(concepts::value_require<Ret>(std::declval<F>()(std::declval<Args>()...)))>
    erased_closure_ref(F& f)
      : _fun_ptr([](void* closure_ptr, Args... args) {
          F& closure = *reinterpret_cast<F*>(closure_ptr);
          return closure(std::forward<Args>(args)...);
      })
      , _closure_ptr(reinterpret_cast<void*>(&f))
    {}

    Ret operator()(Args... args) const { return _fun_ptr(_closure_ptr, std::forward<Args>(args)...); }
};

// void(Args) specialization
template<typename... Args> class erased_closure_ref<void(Args...)>
{
    void (*_fun_ptr)(void*, Args...);
    void* _closure_ptr;

  public:
    erased_closure_ref()                          = delete;
    erased_closure_ref(const erased_closure_ref&) = default;
    erased_closure_ref(erased_closure_ref&)       = default;
    erased_closure_ref& operator=(const erased_closure_ref&) = default;
    erased_closure_ref& operator=(erased_closure_ref&&) = default;

    template<typename F, typename = decltype(std::declval<F>()(std::declval<Args>()...))>
    erased_closure_ref(F& f)
      : _fun_ptr([](void* closure_ptr, Args... args) {
          F& closure = *reinterpret_cast<F*>(closure_ptr);
          closure(std::forward<Args>(args)...);
      })
      , _closure_ptr(reinterpret_cast<void*>(&f))
    {}

    void operator()(Args... args) const { _fun_ptr(_closure_ptr, std::forward<Args>(args)...); }
};

} // namespace gatbl

#endif // GATBL_TYPE_ERASURE
