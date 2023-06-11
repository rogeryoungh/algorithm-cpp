#ifndef ALGO_MATH_POLY_LN
#define ALGO_MATH_POLY_LN

#include "../../base.hpp"
#include "../constant/prepare-inv.hpp"

#include <span>
#include <vector>

template <class ModT, auto poly_div>
auto poly_ln(std::span<const ModT> f, u32 m) {
  std::vector<ModT> x(m);
  std::copy(f.begin(), std::min(f.begin() + m, f.end()), x.begin());
  for (u32 i = 0; i < m; ++i)
    x[i] *= i;
  x = poly_div(x, f, m);
  auto &iv = prepare_inv<ModT>(m);
  for (u32 i = 0; i < m; ++i)
    x[i] *= iv[i];
  return x;
}

#endif // ALGO_MATH_POLY_LN
