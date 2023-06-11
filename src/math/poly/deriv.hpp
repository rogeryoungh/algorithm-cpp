#ifndef ALGO_MATH_POLY_DERIV
#define ALGO_MATH_POLY_DERIV

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <span>
#include <vector>

template <class ModT>
auto poly_deriv(std::span<const ModT> f, u32 m) {
  std::vector<ModT> x(m);
  std::copy(f.begin() + 1, std::min(f.begin() + m + 1, f.end()), x.begin());
  for (u32 i = 0; i < m; ++i)
    x[i] *= i + 1;
  return x;
}

#endif // ALGO_MATH_POLY_DERIV
