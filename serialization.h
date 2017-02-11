#pragma once

#include <vector>
#include <iostream>

template <class T>
std::streamsize binary_size(const std::vector<T>& v)
{
    return v.size() * sizeof(T);
}

template <class T>
void write(std::ostream& s, const T& v)
{
    s.write(reinterpret_cast<const char*>(&v), sizeof(v));
}

template <class T>
void write(std::ostream& s, const std::vector<T>& v)
{
    s.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
    // does not compile with clang due to compiler's bug
    // s.write(reinterpret_cast<const char*>(v.data()), binary_size(v));
}

template <class T>
void read(std::istream& s, T& v)
{
    s.read(reinterpret_cast<char*>(&v), sizeof(v));
}

template <class T>
void read(std::istream& s, std::vector<T>& v)
{
    s.read(reinterpret_cast<char*>(v.data()), binary_size(v));
}

