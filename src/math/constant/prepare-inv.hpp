#ifndef ALGO_MATH_CONSTANT_PREPARE_INV
#define ALGO_MATH_CONSTANT_PREPARE_INV

#include "../../base.hpp"

#include "../../other/modint/modint-concept.hpp"
#include <vector>

template <class ModT>
auto &prepare_inv(u32 m) {
  static std::vector<ModT> iv{1, 1};
  const auto P = ModT::mod();
  while (iv.size() < m) {
    u32 l = iv.size();
    iv.resize(l * 2);
    for (u32 i = l; i < l * 2; ++i)
      iv[i] = iv[P % i] * (P - P / i);
  }
  return iv;
}

#endif // ALGO_MATH_CONSTANT_PREPARE_INV
