#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <cxxabi.h>

#include <typeinfo>
#include <string>
#include <stdexcept>
#include <vector>

/**
 * @brief Return demangled typename of T
 * @see <a href="https://gist.github.com/o11c/3718491">from this gist</a>
 */
template<typename T>
std::string
type_name()
{
    int   status;
    char* result = abi::__cxa_demangle(typeid(T).name(), nullptr, 0, &status);
    switch (status) {
        case 0: break;
        case -1: throw std::bad_alloc();
        case -2: throw std::invalid_argument("bad name argument");
        case -3: throw std::invalid_argument("bad other argument?");
        default: abort();
    }
    try {
        std::string out = result;
        free(result);
        return out;
    } catch (...) {
        free(result);
        throw;
    }
}

static inline uint64_t
splitmix64_stateless(uint64_t index)
{
    uint64_t z = (index + UINT64_C(0x9E3779B97F4A7C15));
    z          = (z ^ (z >> 30u)) * UINT64_C(0xBF58476D1CE4E5B9);
    z          = (z ^ (z >> 27u)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31u);
}

template<typename T = uint64_t>
static inline std::vector<T>
make_randint_vector(size_t n, uint64_t bits, size_t seed = 1)
{
    std::vector<T> res(n);
    const uint64_t mask = (uint64_t(1) << bits) - 1;
    for (size_t i = 0; i < n; i++)
        res[i] = T(splitmix64_stateless(i * seed) & mask);
    return res;
}

#endif // TEST_UTILS_HPP
