#ifndef ALGO_MATH_CONSTANT_PREPARE_INV
#define ALGO_MATH_CONSTANT_PREPARE_INV

#include "../../base.hpp"
#include "../../other/align-alloc.hpp"

template <class ModT>
auto &prepare_inv(u32 m) {
  static AVec<ModT> iv{1, 1};
  const auto P = ModT::mod();
  m = std::bit_ceil(m);
  if (iv.size() < m) {
    u32 n = iv.size();
    iv.resize(m);
    for (u32 i = n; i < m; ++i) {
      iv[i] = iv[P % i] * (P - P / i);
    }
  }
  return iv;
}

#endif // ALGO_MATH_CONSTANT_PREPARE_INV
