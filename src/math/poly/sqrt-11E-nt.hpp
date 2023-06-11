#ifndef ALGO_MATH_POLY_SQRT11ENT
#define ALGO_MATH_POLY_SQRT11ENT

#include "poly-def.hpp"

template <class ModT>
AVec<ModT> poly_sqrt_11E(std::span<const ModT> self, u32 m, const ModT &x0) {
  u32 n = std::bit_ceil(m);
  AVec<ModT> x(n), g(n), ng(n);
  x[0] = x0;
  if (n == 1)
    return x;
  ng[0] = g[0] = x[0].inv();
  x[1] = (self[1] * g[0]).shift2();
  ntt<ModT>({ng.begin(), 2});
  for (u32 t = 2; t < n; t *= 2) {
    AVec<ModT> f(t * 2), nf(t);
    std::copy_n(x.begin(), t, nf.begin());
    ntt<ModT>(nf); // 1E
    std::copy_n(nf.begin(), t, f.begin());
    dot<ModT>(nf, ng);
    intt<ModT>(nf); // 1E
    std::fill_n(nf.begin(), t / 2, 0);
    ntt<ModT>(nf); // 1E
    dot<ModT>(nf, ng);
    intt<ModT>(nf); // 1E
    for (u32 i = t / 2; i < t; ++i)
      g[i] = -nf[i];
    dot<ModT>({f.begin(), t}, f);
    intt<ModT>({f.begin(), t}); // 1E
    for (u32 i = t; i < std::min<u32>(self.size(), t * 2); ++i)
      f[i] = self[i - t] + self[i] - f[i - t];
    std::fill_n(f.begin(), t, 0);
    std::copy_n(g.begin(), t, ng.begin());
    ntt<ModT>(f);                   // 2E
    ntt<ModT>({ng.begin(), t * 2}); // 2E
    dot<ModT>(f, ng);
    intt<ModT>(f); // 2E
    for (u32 i = t; i < t * 2; ++i)
      x[i] = f[i].shift2();
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_SQRT11ENT
