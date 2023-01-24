#ifndef ALGO_MATH_CONSTANT_PREPARE_INV
#define ALGO_MATH_CONSTANT_PREPARE_INV

#include "../../base.hpp"

#include "../../other/modint/modint-concept.hpp"
#include <vector>

template <static_modint_concept ModT>
auto &prepare_inv(u32 m) {
  static std::vector<ModT> iv{1, 1};
  static constexpr auto P = ModT::get_mod();
  while (iv.size() < m) {
    u32 l = iv.size();
    iv.resize(l * 2);
    for (u32 i = l; i < l * 2; ++i)
      iv[i] = iv[P % i] * (P - P / i);
  }
  return iv;
}

#endif // ALGO_MATH_CONSTANT_PREPARE_INV
