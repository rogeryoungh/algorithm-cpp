#ifndef ALGO_H_MATH_AVX2_MONT32_VEC_DOTS
#define ALGO_H_MATH_AVX2_MONT32_VEC_DOTS

#include "../../base.hpp"
#include "../../modular/avx2/mont32x8.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
void dot(ModT *f, const ModT *g, u32 n) {
  if (n <= 16) {
    for (u32 i = 0; i != n; ++i)
      f[i] *= g[i];
  } else {
    using M32x8 = struct M32x8<ModT>;
    auto *fx = reinterpret_cast<M32x8 *>(f);
    auto *gx = reinterpret_cast<const M32x8 *>(g);
    for (u32 i = 0; i != n / 8; ++i)
      fx[i] *= gx[i];
  }
}

template <class ModT>
void div2n(ModT *f, u32 n) {
  ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
  if (n <= 16) {
    for (u32 i = 0; i != n; ++i)
      f[i] *= ivn;
  } else {
    using M32x8 = struct M32x8<ModT>;
    auto *fx = reinterpret_cast<M32x8 *>(f);
    M32x8 ivnx = M32x8::from(ivn);
    for (u32 i = 0; i != n / 8; ++i)
      fx[i] *= ivnx;
  }
}
ALGO_END_NAMESPACE

#endif
