#pragma once

#include "../other/avx512f-helper.hpp"
#include "./montgomery.hpp"

ALGO_BEGIN_NAMESPACE

struct Mont32x16 {
  Mont32 M;
  u32x16 IR, R2, MOD, MOD2, ONE;

  static u32x16 load(const u32 *p) {
    return _mm512_load_si512(reinterpret_cast<const i512 *>(p));
  }

  static void store(u32 *p, u32x16 v) {
    _mm512_store_si512(reinterpret_cast<i512 *>(p), v);
  }

  static u32x16 set1(u32 v) {
    return _mm512_set1_epi32(v);
  }

  Mont32x16(Mont32 mod) : M(mod) {
    IR = set1(mod.IR), R2 = set1(mod.R2);
    MOD = set1(mod.MOD), MOD2 = set1(mod.MOD2);
    ONE = set1(mod.ONE);
  }

  u32x16 norm(u32x16 r) const {
    u32x16 rm = _mm512_sub_epi32(r, MOD);
    return _mm512_min_epu32(r, rm);
  }

  u32x16 add(u32x16 a, u32x16 b) const {
    u32x16 v1 = _mm512_add_epi32(a, b);
    u32x16 v2 = _mm512_sub_epi32(v1, MOD2);
    return _mm512_min_epu32(v1, v2);
  }

  u32x16 sub(u32x16 a, u32x16 b) const {
    u32x16 v1 = _mm512_sub_epi32(a, b);
    u32x16 v2 = _mm512_add_epi32(v1, MOD2);
    return _mm512_min_epu32(v1, v2);
  }

  template <i32 imm>
  u32x16 neg(u32x16 a) const {
    return u32x16_blend<imm>(a, _mm512_sub_epi32(MOD2, a));
  }

  u32x16 reduce(u64x8 x0246, u64x8 x1357) const {
    // (x + u64(u32(x) * IR) * MOD) >> 32;
    auto y0246 = _mm512_mul_epu32(_mm512_mul_epu32(x0246, IR), MOD);
    auto y1357 = _mm512_mul_epu32(_mm512_mul_epu32(x1357, IR), MOD);
    auto z0246 = _mm512_add_epi64(x0246, y0246);
    z0246 = u32x16_permute2301(z0246);
    auto z1357 = _mm512_add_epi64(x1357, y1357);
    return u32x16_blend<0xaaaa>(z0246, z1357);
  }

  u32x16 mul(u32x16 a, u32x16 b) const {
    // return reduce(u64(a) * b);
    u64x8 x0246 = _mm512_mul_epu32(a, b);
    a = u32x16_permute2301(a);
    b = u32x16_permute2301(b);
    u64x8 x1357 = _mm512_mul_epu32(a, b);
    return reduce(x0246, x1357);
  }

  u32x16 trans(u32x16 v) const {
    return mul(v, R2);
  }

  u32x16 get(u32x16 v) const {
    const u32x16 one = set1(1);
    u32x16 v1 = mul(v, one);
    return norm(v1);
  }
};

ALGO_END_NAMESPACE
