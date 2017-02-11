#pragma once

#include <iterator>

#define NDEBUG

#ifndef NDEBUG
#define LOG_DEBUG(expr) std::cerr << expr << '\n';
#else
#define LOG_DEBUG(expr)
#endif

#define LOG_FIELD(field) LOG_DEBUG(#field": " << field)

template <class T>
void print_vector(const std::vector<T>& v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cerr, " "));
    std::cerr << '\n';
}

