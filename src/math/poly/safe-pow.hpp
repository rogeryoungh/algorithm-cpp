#ifndef ALGO_MATH_POLY_SAFE_POW
#define ALGO_MATH_POLY_SAFE_POW

#include "poly-def.hpp"
#include "../constant/prepare-inv.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

template <class ModT, auto poly_pow>
auto poly_safe_pow(std::span<const ModT> f, u64 k, u64 k2, u32 m) {
  auto it = f.begin();
  while (it != f.end() && *it == 0)
    ++it;
  u32 len = it - f.begin();
  if (it == f.end() || (k > 1E9 && len >= 1) || k * len >= m) {
    AVec<ModT> r(m);
    r[0] = k == 0;
    return r;
  }
  AVec<ModT> x(it, f.end());
  ModT f0 = x[0], f0iv = f0.inv(), f0k = f0.pow(k2);

  u32 n = x.size(), i = 0;
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    auto *x8 = reinterpret_cast<X8 *>(x.data());
    u32 nx8 = n / 8;
    X8 f0ivx8 = X8::from(f0iv);
    for (; i < nx8; ++i) {
      x8[i] *= f0ivx8;
    }
  }
#endif
  for (i *= 8; i < n; ++i)
    x[i] *= f0iv;
  x = poly_pow(std::move(x), k, m - k * len);
  n = x.size(), i = 0;
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    auto *x8 = reinterpret_cast<X8 *>(x.data());
    u32 nx8 = n / 8;
    X8 f0kx8 = X8::from(f0k);
    for (; i < nx8; ++i) {
      x8[i] *= f0kx8;
    }
  }
#endif
  for (i *= 8; i < n; ++i)
    x[i] *= f0k;
  x.insert(x.begin(), k * len, 0);
  return x;
}

#endif // ALGO_MATH_POLY_SAFE_POW
