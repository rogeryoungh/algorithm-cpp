#ifndef ALGO_AVX
#define ALGO_AVX

#include "../base.hpp"

// https://judge.yosupo.jp/submission/92714

#pragma GCC target("avx2")

#include <immintrin.h>
#include <array>

namespace simd {

using I256 = __m256i;
using I128x2 = I256;
using I64x4 = I256;
using I32x8 = I256;
using U32x8 = I256;
using U64x4 = I256;
using U128x2 = I256;

namespace i256 {

inline I256 loadu(const I256 *p) {
  return _mm256_loadu_si256(p);
}

template <bool aligned = true>
inline I256 load(const I256 *p) {
  if constexpr (aligned)
    return _mm256_load_si256(p);
  else
    return _mm256_loadu_si256(p);
}

template <bool aligned = true>
inline void store(I256 *p, const I256 &v) {
  if constexpr (aligned)
    _mm256_store_si256(p, v);
  else
    _mm256_storeu_si256(p, v);
}

inline void storeu(I256 *p, const I256 &v) {
  _mm256_storeu_si256(p, v);
}

template <class T>
inline auto to_array(const I256 &v) {
  constexpr u32 sizeT = sizeof(T);
  static_assert(sizeof(I256) % sizeT == 0);
  alignas(32) std::array<T, sizeof(I256) / sizeT> arr;
  _mm256_storeu_si256((I256 *)arr.data(), v);
  return arr;
}

inline I256 bit_and(const I256 &a, const I256 &b) {
  return _mm256_and_si256(a, b);
}

} // namespace i256

namespace i128x2 {

template <i32 imm>
inline I128x2 permute(const I128x2 &a, const I128x2 &b) {
  return _mm256_permute2x128_si256(a, b, imm);
}

template <i32 imm>
inline I128x2 shuffle(const I128x2 &a) {
  return permute<imm>(a, a);
}

} // namespace i128x2

namespace i64x4 {

inline I64x4 add(const I64x4 &a, const I64x4 &b) {
  return _mm256_add_epi64(a, b);
}

} // namespace i64x4

namespace i32x8 {

inline I32x8 from(i32 v) {
  return _mm256_set1_epi32(v);
}

inline I32x8 add(const I32x8 &a, const I32x8 &b) {
  return _mm256_add_epi32(a, b);
}

inline I32x8 sub(const I32x8 &a, const I32x8 &b) {
  return _mm256_sub_epi32(a, b);
}

inline I64x4 mul(const I32x8 &a, const I32x8 &b) {
  return _mm256_mul_epi32(a, b);
}

template <i32 imm>
inline I32x8 shuffle(const I32x8 &a) {
  return _mm256_shuffle_epi32(a, imm);
}

template <i32 imm>
inline I32x8 blend(const I32x8 &a, const I32x8 &b) {
  return _mm256_blend_epi32(a, b, imm);
}

inline I32x8 zero() {
  return _mm256_setzero_si256();
}

inline I32x8 sign(const I32x8 &a) {
  return _mm256_cmpgt_epi32(zero(), a);
}

inline auto mul_0246_1357(const I32x8 &a, const I32x8 &b) {
  auto x0246 = mul(a, b);
  auto x1357 = mul(shuffle<0b11110101>(a), shuffle<0b11110101>(b));
  alignas(32) std::pair<I64x4, I64x4> p = {x0246, x1357};
  return p;
}

inline I32x8 abs(const I32x8 &a) {
  return _mm256_abs_epi32(a);
}

} // namespace i32x8

namespace u32x8 {

inline U64x4 mul(const U32x8 &a, const U32x8 &b) {
  return _mm256_mul_epu32(a, b);
}

inline auto mul_0246_1357(const U32x8 &a, const U32x8 &b) {
  auto x0246 = mul(a, b);
  auto x1357 = mul(i32x8::shuffle<0b11110101>(a), i32x8::shuffle<0b11110101>(b));
  alignas(32) std::pair<U64x4, U64x4> p = {x0246, x1357};
  return p;
}

} // namespace u32x8

} // namespace simd

#endif
