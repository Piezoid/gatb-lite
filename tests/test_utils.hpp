#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <cxxabi.h>

#include <typeinfo>
#include <string>
#include <stdexcept>

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

#endif // TEST_UTILS_HPP
