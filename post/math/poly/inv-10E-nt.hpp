#ifndef ALGO_MATH_POLY_INV10ENT
#define ALGO_MATH_POLY_INV10ENT

#include "poly-def.hpp"

template <class ModT>
auto poly_inv_10E(std::span<const ModT> f, u32 m) {
  u32 n = std::bit_ceil(m);
  AVec<ModT> x(n);
  x[0] = f[0].inv();
  for (u32 t = 1; t < n; t *= 2) {
    AVec<ModT> f2(t * 2), nx(t * 2);
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
    vectorization_2<ModT, true>(t, x.data() + t, f2.data() + t, []<class T>(T &xi, T f2) {
      xi = -f2;
    });
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INV10ENT
