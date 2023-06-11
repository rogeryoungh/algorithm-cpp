#ifndef ALGO_MATH_POLY_INVSQRT_12ENT
#define ALGO_MATH_POLY_INVSQRT_12ENT

#include "poly-def.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

template <class ModT>
auto poly_invsqrt_12E(std::span<const ModT> self, u32 m, const ModT &x0) {
  u32 n = std::bit_ceil(m);
  AVec<ModT> x(n * 2);
  x[0] = x0;
  ModT ivn2 = -ModT(2).inv();
  for (u32 t = 1; t < n; t *= 2) {
    AVec<ModT> u(t * 4), s(t * 4);
    std::copy(self.begin(), std::min(self.begin() + t * 2, self.end()), u.begin());
    std::copy(x.begin(), x.begin() + t, s.begin());
    ntt<ModT>(u); // 4E
    ntt<ModT>(s); // 4E
#ifndef ALGO_DISABLE_SIMD_AVX2
    if (montgomery_modint_concept<ModT> && t * 4 > 16) {
      using X8 = simd::M32x8<ModT>;
      auto *ux8 = reinterpret_cast<X8 *>(u.data());
      auto *sx8 = reinterpret_cast<X8 *>(s.data());
      X8 ivn2x8 = X8::from(ivn2);
      u32 nx8 = t * 4 / 8;
      for (u32 i = 0; i < nx8; ++i) {
        ux8[i] *= sx8[i] * sx8[i] * sx8[i] * ivn2x8;
      }
    } else {
#endif
      for (u32 i = 0; i < t * 4; ++i) {
        u[i] *= s[i] * s[i] * s[i] * ivn2;
      }
#ifndef ALGO_DISABLE_SIMD_AVX2
    }
#endif
    intt<ModT>(u); // 4E
    std::copy_n(u.begin() + t, t, x.begin() + t);
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INVSQRT_12ENT
