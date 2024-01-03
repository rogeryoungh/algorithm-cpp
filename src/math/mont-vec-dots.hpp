#ifndef ALGO_H_MATH_MONT32_VEC_DOTS
#define ALGO_H_MATH_MONT32_VEC_DOTS

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
void dot(ModT *f, const ModT *g, u32 n) {
  for (u32 i = 0; i != n; ++i)
    f[i] *= g[i];
}

template <class ModT>
void div2n(ModT *f, u32 n) {
  ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
  for (u32 i = 0; i != n; ++i)
    f[i] *= ivn;
}
ALGO_END_NAMESPACE

#endif
