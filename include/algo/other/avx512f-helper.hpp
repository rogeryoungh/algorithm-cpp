#pragma once

#include "./use-avx512f.hpp"

ALGO_BEGIN_NAMESPACE

template <_MM_PERM_ENUM imm>
inline u32x16 u32x16_shuffle(u32x16 a) {
  return _mm512_shuffle_epi32(a, imm);
}

template <__mmask16 imm>
inline u32x16 u32x16_blend(u32x16 a, u32x16 b) {
  return _mm512_mask_blend_epi32(imm, a, b);
}

inline u32x16 u32x16_permute2301(u32x16 a) { // 1, 0.5
  return u32x16_shuffle<_MM_PERM_DDBB>(a);
}

ALGO_END_NAMESPACE
