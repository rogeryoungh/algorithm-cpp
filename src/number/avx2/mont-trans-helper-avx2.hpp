#ifndef ALGO_H_MATH_MONT32_TRANS_HELPER
#define ALGO_H_MATH_MONT32_TRANS_HELPER

#include "../../base.hpp"
#include "./mont32x8.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
void mont32_trans_avx2(ModT *out, const u32 *in, u32 n) {
  if (n < 8) {
    for (u32 i = 0; i != n; ++i)
      out[i] = in[i];
  } else {
    M32x8<ModT>::trans(out, in, n);
  }
}

template <class ModT>
void mont32_get_avx2(u32 *out, const ModT *in, u32 n) {
  if (n < 8) {
    for (u32 i = 0; i != n; ++i)
      out[i] = in[i].get();
  } else {
    M32x8<ModT>::get(out, in, n);
  }
}

ALGO_END_NAMESPACE

#endif
