#ifndef ALGO_H_MATH_MONT32_VEC_DOTS
#define ALGO_H_MATH_MONT32_VEC_DOTS

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
void mont_dot(ModT *f, const ModT *g, u32 n) {
  for (u32 i = 0; i != n; ++i)
    f[i] *= g[i];
}

template <class ModT>
void mont_dot1(ModT *f, u32 n, const ModT g0) {
  for (u32 i = 0; i != n; ++i)
    f[i] *= g0;
}

ALGO_END_NAMESPACE

#endif
