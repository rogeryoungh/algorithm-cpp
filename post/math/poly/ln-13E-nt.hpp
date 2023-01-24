#ifndef ALGO_MATH_POLY_LN13ENT
#define ALGO_MATH_POLY_LN13ENT

#include "../../base.hpp"
#include "div-13E-nt.hpp"
#include "../constant/prepare-inv.hpp"

#include <span>
#include <vector>

template <static_modint_concept ModT>
auto poly_ln_13E(std::span<const ModT> f, u32 m) {
  std::vector<ModT> x(m);
  std::copy(f.begin(), std::min(f.begin() + m, f.end()), x.begin());
  for (u32 i = 0; i < m; ++i)
    x[i] *= i;
  x = poly_div_13E<ModT>(x, f, m);
  auto &iv = prepare_inv<ModT>(m);
  for (u32 i = 0; i < m; ++i)
    x[i] *= iv[i];
  return x;
}

#endif // ALGO_MATH_POLY_LN13ENT
