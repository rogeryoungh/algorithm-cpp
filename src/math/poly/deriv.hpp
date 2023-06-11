#ifndef ALGO_MATH_POLY_DERIV
#define ALGO_MATH_POLY_DERIV

#include "poly-def.hpp"

template <class ModT>
auto poly_deriv(std::span<const ModT> f, u32 m) {
  AVec<ModT> x(m);
  std::copy(f.begin() + 1, std::min(f.begin() + m + 1, f.end()), x.begin());
  for (u32 i = 0; i < m; ++i)
    x[i] *= i + 1;
  return x;
}

#endif // ALGO_MATH_POLY_DERIV
