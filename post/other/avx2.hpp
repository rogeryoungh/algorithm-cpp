#ifndef ALGO_SIMD_AVX2
#define ALGO_SIMD_AVX2

#include "../base.hpp"

// https://judge.yosupo.jp/submission/92714

#pragma GCC target("avx2")

#include <immintrin.h>
#include <array>

namespace simd {

using i256 = __m256i;
using i256u = __m256i_u;

template <bool aligned>
using i256a = std::conditional_t<aligned, simd::i256, simd::i256u>;

using i32x8 = __m256i;
using i64x4 = __m256i;
using u32x8 = __m256i;
using u64x4 = __m256i;

template <class T>
inline auto i256_toarray(const i256 &v) {
  constexpr u32 sizeT = sizeof(T);
  static_assert(sizeof(i256) % sizeT == 0);
  alignas(32) std::array<T, sizeof(i256) / sizeT> arr;
  _mm256_storeu_si256((i256 *)arr.data(), v);
  return arr;
}

inline i32x8 i32x8_set1(i32 v) {
  return _mm256_set1_epi32(v);
}

inline i32x8 i32x8_sign(i32x8 a) { // 2, 0.33
  // v = (v >> 31) ? -1 : 0;
  return _mm256_cmpgt_epi32(_mm256_setzero_si256(), a);
}

inline u64x4 u32x8_mul0246(u32x8 a, u32x8 b) { // 5, 0.5
  return _mm256_mul_epu32(a, b);
}

inline i32x8 i32x8_swap_lohi(i32x8 a) { // 1, 0.5
  return _mm256_shuffle_epi32(a, 0b11110101);
}

inline u64x4 u32x8_mul1357(u32x8 a, u32x8 b) { // 7, 0.5
  a = i32x8_swap_lohi(a);
  b = i32x8_swap_lohi(b);
  return u32x8_mul0246(a, b);
}

} // namespace simd

#endif // ALGO_SIMD_AVX2
