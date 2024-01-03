#ifndef ALGO_H_MATH_AVX512F_MONT32_VEC_DOTS
#define ALGO_H_MATH_AVX512F_MONT32_VEC_DOTS

#include "../../base.hpp"
#include "../../modular/avx512f/mont32x16.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
void dot(ModT *f, const ModT *g, u32 n) {
  using M32x16 = struct M32x16<ModT>;
  if (n <= 16) {
    for (u32 i = 0; i != n; ++i)
      f[i] *= g[i];
  } else {
    auto *fx = reinterpret_cast<M32x16 *>(f);
    auto *gx = reinterpret_cast<const M32x16 *>(g);
    for (u32 i = 0; i != n / 16; ++i)
      fx[i] *= gx[i];
  }
}

template <class ModT>
void div2n(ModT *f, u32 n) {
  using M32x16 = struct M32x16<ModT>;
  ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
  if (n <= 16) {
    for (u32 i = 0; i != n; ++i)
      f[i] *= ivn;
  } else {
    auto *fx = reinterpret_cast<M32x16 *>(f);
    M32x16 invx = M32x16::from(ivn);
    for (u32 i = 0; i != n / 16; ++i)
      fx[i] *= invx;
  }
}
ALGO_END_NAMESPACE

#endif
