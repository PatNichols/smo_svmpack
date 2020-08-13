#ifndef PUTILS_HPP
#define PUTILS_HPP
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstddef>

namespace putils {


inline void ToUpper(std::string& s) {
    for (size_t k=0;k<s.size();++k) s[k]=toupper(s[k]);
}


template < typename T > 
inline T StringToType(const std::string& str) {
    try {
        std::istringstream is(str);
        T x;
        is >> x;
        return x;
    } catch (std::exception& e) {
        std::cerr << " StringToType cannot convert " << str << " to type requested\n";    
        exit(EXIT_FAILURE);
    }
}

template <> inline int StringToType<int>(const std::string& str) { return std::stoi(str);}
template <> inline unsigned int StringToType<unsigned int>(const std::string& str) { return static_cast<unsigned int>(std::stoi(str));}
template <> inline long StringToType<long>(const std::string& str) { return std::stol(str);}
template <> inline unsigned long StringToType<unsigned long>(const std::string& str) { return std::stoul(str);}
template <> inline long long int StringToType<long long int>(const std::string& str) { return std::stoll(str);}
template <> inline unsigned long long int StringToType<unsigned long long int>(const std::string& str) 
{ return std::stoull(str);}

template <> inline float StringToType<float>(const std::string& str) { return std::stof(str);}
template <> inline double StringToType<double>(const std::string& str) { return std::stod(str);}
template <> inline long double StringToType<long double>(const std::string& str) { return std::stold(str);}

template <> inline const char* StringToType< const char * >(const std::string& str) { return str.c_str();}

template <> inline bool StringToType<bool>(const std::string& str) {
    if (str.size()==0) {
        std::cerr << "empty string passed to StringToType\n";
        exit(EXIT_FAILURE); 
    }
    char ch = str[0];
    if (ch=='N' || ch=='n' || ch=='0' || ch=='F' || ch=='f') {
        return false;
    }    
    return true;
}

template <> inline short StringToType<short>(const std::string& str) {
    int x = stoi(str);
    return (short)x;
}

template <> inline unsigned short StringToType<unsigned short>(const std::string& str) {
    int x = stoi(str);
    return (unsigned short)x;
}

template < typename T > inline std::string&& TypeToString( const T& x) {
    std::ostringstream out;
    out << x;    
    return std::move(out.str());
}

} // end namespace
#endif