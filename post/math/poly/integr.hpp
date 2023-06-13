#ifndef ALGO_MATH_POLY_INTEGR
#define ALGO_MATH_POLY_INTEGR

#include "poly-def.hpp"
#include "../constant/prepare-inv.hpp"

template <class ModT>
auto poly_integr(std::span<const ModT> f, u32 m, u32 C) {
  AVec<ModT> x(m);
  std::copy(f.begin(), std::min(f.begin() + m - 1, f.end()), x.begin() + 1);
  const auto &iv = prepare_inv<ModT>(m);
  dot<ModT>(x, iv);
  x[0] = C;
  return x;
}

#endif // ALGO_MATH_POLY_DIV13ENT
