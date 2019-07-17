#ifndef GATBL_STRING_HPP
#define GATBL_STRING_HPP

#include <string>
#include <ostream>

namespace gatbl {

inline bool
hasEnding(std::string const& fullString, std::string const& ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

template<typename T>
auto
operator<<(std::ostream& out, T x) -> decltype(out << to_string(x))
{
    return out << to_string(x);
}

} // namespace gatbl

#endif // GATBL_STRING_HPP
