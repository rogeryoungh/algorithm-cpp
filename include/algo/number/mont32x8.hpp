#pragma once

#include "../other/avx2-helper.hpp"
#include "./montgomery.hpp"

ALGO_BEGIN_NAMESPACE

struct Mont32x8 {
  Mont32 M;
  u32x8 IR, R2, MOD, MOD2, ONE;

  static u32x8 loadu(const u32 *p) {
    return _mm256_loadu_si256(reinterpret_cast<const i256 *>(p));
  }

  static void storeu(u32 *p, u32x8 v) {
    _mm256_storeu_si256(reinterpret_cast<i256 *>(p), v);
  }

  static u32x8 set1(u32 v) {
    return _mm256_set1_epi32(v);
  }

  Mont32x8(Mont32 mod) : M(mod) {
    IR = set1(mod.IR), R2 = set1(mod.R2);
    MOD = set1(mod.MOD), MOD2 = set1(mod.MOD2);
    ONE = set1(mod.ONE);
  }

  u32x8 norm(u32x8 r) const {
    u32x8 rm = _mm256_sub_epi32(r, MOD);
    return _mm256_min_epu32(r, rm);
  }

  u32x8 add(u32x8 a, u32x8 b) const {
    u32x8 v1 = _mm256_add_epi32(a, b);
    u32x8 v2 = _mm256_sub_epi32(v1, MOD2);
    return _mm256_min_epu32(v1, v2);
  }

  u32x8 sub(u32x8 a, u32x8 b) const {
    u32x8 v1 = _mm256_sub_epi32(a, b);
    u32x8 v2 = _mm256_add_epi32(v1, MOD2);
    return _mm256_min_epu32(v1, v2);
  }

  template <i32 imm>
  u32x8 neg(u32x8 a) const {
    return u32x8_blend<imm>(a, _mm256_sub_epi32(MOD2, a));
  }

  u32x8 reduce(u64x4 x0246, u64x4 x1357) const {
    // (x + u64(u32(x) * IR) * MOD) >> 32;
    auto y0246 = _mm256_mul_epu32(_mm256_mul_epu32(x0246, IR), MOD);
    auto y1357 = _mm256_mul_epu32(_mm256_mul_epu32(x1357, IR), MOD);
    auto z0246 = _mm256_add_epi64(x0246, y0246);
    z0246 = u32x8_permute2301(z0246);
    auto z1357 = _mm256_add_epi64(x1357, y1357);
    return u32x8_blend<0xaa>(z0246, z1357);
  }

  u32x8 mul(u32x8 a, u32x8 b) const {
    // return reduce(u64(a) * b);
    u64x4 x0246 = _mm256_mul_epu32(a, b);
    a = u32x8_permute2301(a);
    b = u32x8_permute2301(b);
    u64x4 x1357 = _mm256_mul_epu32(a, b);
    return reduce(x0246, x1357);
  }

  u32x8 trans(u32x8 v) const {
    return mul(v, R2);
  }

  u32x8 get(u32x8 v) const {
    const u32x8 one = set1(1);
    u32x8 v1 = mul(v, one);
    return norm(v1);
  }
};

ALGO_END_NAMESPACE
