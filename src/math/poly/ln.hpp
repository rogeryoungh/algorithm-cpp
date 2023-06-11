#ifndef ALGO_MATH_POLY_LN
#define ALGO_MATH_POLY_LN

#include "poly-def.hpp"
#include "../constant/prepare-inv.hpp"

template <class ModT, auto poly_div>
auto poly_ln(std::span<const ModT> f, u32 m) {
  AVec<ModT> x(m);
  std::copy(f.begin(), std::min(f.begin() + m, f.end()), x.begin());
  for (u32 i = 0; i < m; ++i)
    x[i] *= i;
  x = poly_div(x, f, m);
  const auto &iv = prepare_inv<ModT>(m);
  dot<ModT>(x, iv);
  return x;
}

#endif // ALGO_MATH_POLY_LN
