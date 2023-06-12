#ifndef ALGO_MATH_POLY_INV10ENT
#define ALGO_MATH_POLY_INV10ENT

#include "poly-def.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

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
#ifndef ALGO_DISABLE_SIMD_AVX2
    if (montgomery_modint_concept<ModT> && t * 2 > 16) {
      using X8 = simd::M32x8<ModT>;
      auto *x8 = reinterpret_cast<X8 *>(x.data());
      auto *f2x8 = reinterpret_cast<X8 *>(f2.data());
      u32 nx8 = t * 2 / 8;
      for (u32 i = nx8 / 2; i < nx8; ++i) {
        x8[i] = -f2x8[i];
      }
    } else {
#endif
      for (u32 i = t; i < t * 2; ++i) {
        x[i] = -f2[i];
      }
#ifndef ALGO_DISABLE_SIMD_AVX2
    }
#endif
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INV10ENT
