#ifndef ALGO_MATH_POLY_INV10ENT
#define ALGO_MATH_POLY_INV10ENT

#include "../../base.hpp"
#include "ntt.hpp"
#include "vec-dots.hpp"

#include <algorithm>
#include <vector>

template <static_modint_concept ModT>
auto poly_inv_10E(std::span<const ModT> f, u32 m) {
  u32 n = std::bit_ceil(m);
  std::vector<ModT> x(n);
  x[0] = f[0].inv();
  for (u32 t = 1; t < n; t *= 2) {
    std::vector<ModT> f2(t * 2), nx(t * 2);
    std::copy(f.begin(), std::min(f.begin() + t * 2, f.end()), f2.begin());
    std::copy_n(x.begin(), t * 2, nx.begin());
    ntt<ModT>(f2); // 2E
    ntt<ModT>(nx); // 2E
    dot<ModT>(f2, nx);
    intt<ModT>(f2); // 2E
    std::fill_n(f2.begin(), t, 0);
    ntt<ModT>(f2); // 2E
    dot<ModT>(f2, nx);
    intt<ModT>(f2); // 2E
    for (u32 i = t; i < t * 2; ++i) {
      x[i] = -f2[i];
    }
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INV10ENT
