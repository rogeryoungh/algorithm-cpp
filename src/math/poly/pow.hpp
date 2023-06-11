#ifndef ALGO_MATH_POLY_POW
#define ALGO_MATH_POLY_POW

#include "poly-def.hpp"
#include "../constant/prepare-inv.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

template <class ModT, auto poly_ln, auto poly_exp>
auto poly_pow(std::span<const ModT> f, u64 k, u32 m) {
  if (k == 0) {
    AVec<ModT> x(m);
    x[0] = 1;
    return x;
  } else {
    ModT mk = ModT::safe(k);
    auto x = poly_ln(f, m);
    u32 n = x.size(), i = 0;
#ifndef ALGO_DISABLE_SIMD_AVX2
    if (montgomery_modint_concept<ModT> && n > 16) {
      using X8 = simd::M32x8<ModT>;
      auto *x8 = reinterpret_cast<X8 *>(x.data());
      u32 nx8 = n / 8;
      X8 mkx8 = X8::from(mk);
      for (; i < nx8; ++i) {
        x8[i] *= mkx8;
      }
    }
#endif
    for (i *= 8; i < n; ++i)
      x[i] *= mk;
    return poly_exp(std::move(x), m);
  }
}

#endif // ALGO_MATH_POLY_LN
