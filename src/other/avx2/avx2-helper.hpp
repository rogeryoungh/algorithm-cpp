#ifndef ALGO_SIMD_AVX2
#define ALGO_SIMD_AVX2

#include "./use-avx2.hpp"

ALGO_BEGIN_NAMESPACE

inline i256 i256_load(void *p) {
  return _mm256_load_si256(reinterpret_cast<i256 *>(p));
}

template <int imm>
inline u32x8 u32x8_shuffle(u32x8 a) {
  return _mm256_shuffle_epi32(a, imm);
}

template <int imm>
inline u32x8 u32x8_blend(u32x8 a, u32x8 b) {
  return _mm256_blend_epi32(a, b, imm);
}

inline u64x4 u32x8_mul0246(u32x8 a, u32x8 b) { // 5, 0.5
  return _mm256_mul_epu32(a, b);
}

inline i32x8 i32x8_swap_lohi(i32x8 a) { // 1, 0.5
  return u32x8_shuffle<0xf5>(a);
}

inline u64x4 u32x8_mul1357(u32x8 a, u32x8 b) { // 7, 0.5
  a = i32x8_swap_lohi(a);
  b = i32x8_swap_lohi(b);
  return u32x8_mul0246(a, b);
}

ALGO_END_NAMESPACE

#endif
