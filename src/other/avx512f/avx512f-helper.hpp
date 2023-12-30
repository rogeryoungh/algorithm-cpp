#ifndef ALGO_SIMD_AVX512F
#define ALGO_SIMD_AVX512F

#include "./use-avx512f.hpp"

ALGO_BEGIN_NAMESPACE

inline i512 i512_load(void *p) {
  return _mm512_load_si512(reinterpret_cast<i512 *>(p));
}

template <_MM_PERM_ENUM imm>
inline u32x16 u32x16_shuffle(u32x16 a) {
  return _mm512_shuffle_epi32(a, imm);
}

template <__mmask16 imm>
inline u32x16 u32x16_blend(u32x16 a, u32x16 b) {
  return _mm512_mask_blend_epi32(imm, a, b);
}

inline u64x8 u32x16_mul0246(u32x16 a, u32x16 b) { // 5, 0.5
  return _mm512_mul_epu32(a, b);
}

inline i32x16 i32x16_swap_lohi(i32x16 a) { // 1, 0.5
  return u32x16_shuffle<_MM_PERM_DDBB>(a);
}

inline u64x8 u32x16_mul1357(u32x16 a, u32x16 b) { // 7, 0.5
  a = i32x16_swap_lohi(a);
  b = i32x16_swap_lohi(b);
  return u32x16_mul0246(a, b);
}

ALGO_END_NAMESPACE

#endif
