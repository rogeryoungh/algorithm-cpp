#ifndef ALGO_MATH_POLY_INTEGR
#define ALGO_MATH_POLY_INTEGR

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"
#include "../constant/prepare-inv.hpp"

#include <span>
#include <vector>

template <static_modint_concept ModT>
auto poly_integr(std::span<const ModT> f, u32 m, u32 C) {
  std::vector<ModT> x(m);
  std::copy(f.begin(), std::min(f.begin() + m - 1, f.end()), x.begin() + 1);
  auto &iv = prepare_inv<ModT>(m);
  for (u32 i = 1; i < m; ++i)
    x[i] *= iv[i];
  x[0] = C;
  return x;
}

#endif // ALGO_MATH_POLY_DIV13ENT
