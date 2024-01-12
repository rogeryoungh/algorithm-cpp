#ifndef ALGO_H_MATH_MONT32_TRANS_HELPER
#define ALGO_H_MATH_MONT32_TRANS_HELPER

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
void mont_trans(ModT *out, const u32 *in, u32 n) {
  for (u32 i = 0; i != n; ++i)
    out[i] = in[i];
}

template <class ModT>
void mont_get(u32 *out, const ModT *in, u32 n) {
  for (u32 i = 0; i != n; ++i)
    out[i] = in[i].get();
}

ALGO_END_NAMESPACE

#endif
