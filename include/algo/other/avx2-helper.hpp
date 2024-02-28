#pragma once

#include "./use-avx2.hpp"

ALGO_BEGIN_NAMESPACE

template <int imm>
inline u32x8 u32x8_shuffle(u32x8 a) {
  return _mm256_shuffle_epi32(a, imm);
}

template <int imm>
inline u32x8 u32x8_blend(u32x8 a, u32x8 b) {
  return _mm256_blend_epi32(a, b, imm);
}

inline i32x8 i32x8_permute2301(i32x8 a) { // 1, 0.5
  return u32x8_shuffle<0xf5>(a);
}

// https://stackoverflow.com/questions/37296289/fastest-way-to-multiply-an-array-of-int64-t
u128x2 u64x4_mul0246(u64x4 a, u64x4 b) {
  u64x4 b_swap = _mm256_shuffle_epi32(b, _MM_SHUFFLE(2, 3, 0, 1));
  u64x4 crossprod = _mm256_mullo_epi32(a, b_swap);

  u64x4 prodlh = _mm256_slli_epi64(crossprod, 32);
  u64x4 prodhl = _mm256_and_si256(crossprod, _mm256_set1_epi64x(0xFFFFFFFF00000000));
  u64x4 sumcross = _mm256_add_epi32(prodlh, prodhl);

  u64x4 prodll = _mm256_mul_epu32(a, b);
  u64x4 prod = _mm256_add_epi32(prodll, sumcross);
  return prod;
}

ALGO_END_NAMESPACE
