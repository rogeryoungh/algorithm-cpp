#ifndef ALGO_H_MODULAR_MONT32X8
#define ALGO_H_MODULAR_MONT32X8

#include "../../base.hpp"
#include "../../other/avx2/use-avx2.hpp"
#include "../../other/avx2/avx2-helper.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct M32x8 {
  u32x8 v;

  M32x8() : v() {}
  M32x8(const i256 &o) : v(o) {}

  operator u32x8() const {
    return v;
  }

  static u32x8 from(ModT v) {
    return _mm256_set1_epi32(v.raw());
  }

  M32x8 &operator+=(M32x8 r) {
    u32x8 v1 = _mm256_add_epi32(v, r.v);
    u32x8 v2 = _mm256_sub_epi32(v1, mod2x8());
    v = _mm256_min_epu32(v1, v2);
    return *this;
  }

  M32x8 &operator-=(M32x8 r) {
    u32x8 v1 = _mm256_sub_epi32(v, r.v);
    u32x8 v2 = _mm256_add_epi32(v1, mod2x8());
    v = _mm256_min_epu32(v1, v2);
    return *this;
  }

  friend M32x8 operator+(const M32x8 &lhs, const M32x8 &rhs) {
    return M32x8(lhs) += rhs;
  }

  friend M32x8 operator-(const M32x8 &lhs, const M32x8 &rhs) {
    return M32x8(lhs) -= rhs;
  }

  friend M32x8 operator*(const M32x8 &lhs, const M32x8 &rhs) {
    return M32x8(lhs) *= rhs;
  }

  template <i32 imm>
  M32x8  neg() const {
    return u32x8_blend<imm>(v, _mm256_sub_epi32(mod2x8(), v));
  }

  static u32x8 reduce(u64x4 x0246, u64x4 x1357) {
    // (x + u64(u32(x) * IR) * MOD) >> 32;
    auto y0246 = u32x8_mul0246(u32x8_mul0246(x0246, irx8()), modx8());
    auto y1357 = u32x8_mul0246(u32x8_mul0246(x1357, irx8()), modx8());
    auto z0246 = _mm256_add_epi64(x0246, y0246);
    z0246 = i32x8_permute2301(z0246);
    auto z1357 = _mm256_add_epi64(x1357, y1357);
    return u32x8_blend<0xaa>(z0246, z1357);
  }

  static u32x8 mul_reduce(u32x8 a, u32x8 b) {
    // return reduce(u64(a) * b);
    u64x4 x0246 = u32x8_mul0246(a, b);
    u64x4 x1357 = u32x8_mul1357(a, b);
    return reduce(x0246, x1357);
  }

  static u32x8 trans(u32x8 v) {
    return mul_reduce(v, r2x8());
  }

  static u32x8 get(u32x8 v) {
    u32x8 v1 = mul_reduce(v, _mm256_set1_epi32(1));
    u32x8 v2 = _mm256_sub_epi32(v1, modx8());
    return _mm256_min_epu32(v1, v2);
  }

  M32x8 &operator*=(const M32x8 &rhs) {
    v = mul_reduce(v, rhs.v);
    return *this;
  }

  static u32x8 div2(u32x8 v) {
    // (v & 1) ? v + mod : v
    u32x8 one = _mm256_set1_epi32(1);
    u32x8 and1 = _mm256_and_si256(one, v);
    and1 = _mm256_sub_epi32(and1, one);
    u32x8 smodx8 = _mm256_andnot_si256(and1, modx8());
    u32x8 ans = _mm256_add_epi32(v, smodx8);
    return _mm256_srli_epi32(ans, 1);
  }

  static u32x8 irx8() {
    return _mm256_set1_epi32(ModT::IR);
  }
  static u32x8 r2x8() {
    return _mm256_set1_epi32(ModT::R2);
  }
  static u32x8 mod2x8() {
    return _mm256_set1_epi32(ModT::MOD2);
  }
  static u32x8 modx8() {
    return _mm256_set1_epi32(ModT::MOD);
  }
};

ALGO_END_NAMESPACE

#endif
