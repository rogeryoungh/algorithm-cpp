#ifndef ALGO_MATH_POLY_EXP17ENT
#define ALGO_MATH_POLY_EXP17ENT

#include "../../base.hpp"
#include "ntt.hpp"
#include "../constant/prepare-inv.hpp"

#include <algorithm>
#include <vector>
#include <iostream>

template <class ModT>
std::vector<ModT> poly_exp_17E(std::span<const ModT> self, u32 m) {
  u32 n = std::bit_ceil(m);
  auto &iv = prepare_inv<ModT>(n);
  std::vector<ModT> f(n), g(n);
  f[0] = g[0] = 1;
  if (m == 1)
    return f;
  f[1] = self[1];
  auto nf = f, ng = g;
  ntt<ModT>({ng.begin(), 2});

  for (u32 t = 2; t < n; t *= 2) {
    std::copy_n(f.begin(), t, nf.begin());
    ntt<ModT>(std::span(nf.begin(), t * 2)); // 2E

    std::vector<ModT> q(t * 2);
    std::span qt{q.begin(), t};
    std::copy_n(nf.begin(), t * 2, q.begin());

    dot<ModT>(qt, ng);
    intt<ModT>(qt); // 1E
    std::fill_n(q.begin(), t / 2, 0);
    ntt<ModT>(qt); // 1E
    dot<ModT>(qt, ng);
    intt<ModT>(qt); // 1E
    for (u32 i = t / 2; i < t; ++i)
      g[i] = -q[i];
    std::copy_n(g.begin(), t, ng.begin());

    std::fill_n(q.begin() + t, t, 0);
    for (u32 i = 0; i < t; ++i)
      q[i] = self[i] * i;
    ntt<ModT>(qt); // 1E
    dot<ModT>(qt, nf);
    intt<ModT>(qt); // 1E
    for (u32 i = 0; i < t; ++i)
      q[i] -= f[i] * i;
    ntt<ModT>(q);                   // 2E
    ntt<ModT>({ng.begin(), t * 2}); // 2E
    dot<ModT>(q, ng);
    intt<ModT>(q); // 2E
    for (u32 i = t; i < t * 2; ++i)
      q[i] = q[i - t] * iv[i];
    for (u32 i = t; i < std::min<u32>(t * 2, self.size()); ++i)
      q[i] += self[i];
    std::fill_n(q.begin(), t, 0);
    ntt<ModT>(q); // 2E
    dot<ModT>(q, nf);
    intt<ModT>(q); // 2E
    std::copy_n(q.begin() + t, t, f.begin() + t);
  }
  return f.resize(m), f;
}

#endif // ALGO_MATH_POLY_DIV13ENT
