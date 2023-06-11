#ifndef ALGO_MATH_POLY_INV12ENT
#define ALGO_MATH_POLY_INV12ENT

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"
#include "ntt.hpp"

#include <algorithm>
#include <vector>

template <class ModT>
auto poly_inv_12E(std::span<const ModT> f, u32 m) {
  u32 n = std::bit_ceil(m);
  std::vector<ModT> x(n * 2);
  x[0] = f[0].inv();
  for (u32 t = 1; t < n; t *= 2) {
    std::span xt{x.begin(), x.begin() + t * 4};
    std::vector<ModT> u(t * 4);
    std::copy(f.begin(), std::min(f.begin() + t * 2, f.end()), u.begin());
    ntt<ModT>(u);  // 4E
    ntt<ModT>(xt); // 4E
    for (u32 i = 0; i < t * 4; ++i) {
      x[i] = (ModT(2) + -u[i] * x[i]) * x[i];
    }
    intt<ModT>(xt); // 4E
    std::fill_n(x.begin() + t * 2, t * 2, 0);
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INV12ENT
