#ifndef GATBL_NOP_FUNCTOR_HPP
#define GATBL_NOP_FUNCTOR_HPP

namespace gatbl {

struct nop_functor
{
    nop_functor() noexcept                   = default;
    nop_functor(nop_functor&&) noexcept      = default;
    nop_functor(const nop_functor&) noexcept = default;

    template<typename... Args> void operator()(Args...) const noexcept {}
};

}

#endif // GATBL_NOP_FUNCTOR_HPP
