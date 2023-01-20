#ifndef ALGO_DEBUG_HPP
#define ALGO_DEBUG_HPP

#include "../base.hpp"

#include <iostream>
#include <vector>

#define dbg(x) debug(__FILE__, __LINE__, #x, x)

void debug(const std::string &file, i32 line, const std::string &name, const auto &value) {
  std::cerr << file << ":" << line << " | " << name << " = " << value << std::endl;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  os << "[ ";
  for (const T &vi : v) {
    os << vi << ", ";
  }
  return os << "]";
}

#endif // ALGO_DEBUG_HPP
