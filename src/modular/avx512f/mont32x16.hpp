#ifndef ALGO_H_MODULAR_MONT32X16
#define ALGO_H_MODULAR_MONT32X16

#include "../../base.hpp"
#include "../../other/avx512f/use-avx512f.hpp"
#include "../../other/avx512f/avx512f-helper.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct M32x16 {
  u32x16 v;

  M32x16() : v() {}
  M32x16(const i512 &o) : v(o) {}

  operator u32x16() const {
    return v;
  }

  static u32x16 from(ModT v) {
    return _mm512_set1_epi32(v.raw());
  }

  M32x16 &operator+=(M32x16 r) {
    u32x16 v1 = _mm512_add_epi32(v, r.v);
    u32x16 v2 = _mm512_sub_epi32(v1, mod2x16());
    v = _mm512_min_epu32(v1, v2);
    return *this;
  }

  M32x16 &operator-=(M32x16 r) {
    u32x16 v1 = _mm512_sub_epi32(v, r.v);
    u32x16 v2 = _mm512_add_epi32(v1, mod2x16());
    v = _mm512_min_epu32(v1, v2);
    return *this;
  }

  friend M32x16 operator+(const M32x16 &lhs, const M32x16 &rhs) {
    return M32x16(lhs) += rhs;
  }

  friend M32x16 operator-(const M32x16 &lhs, const M32x16 &rhs) {
    return M32x16(lhs) -= rhs;
  }

  friend M32x16 operator*(const M32x16 &lhs, const M32x16 &rhs) {
    return M32x16(lhs) *= rhs;
  }

  template <i32 imm>
  M32x16 neg() const {
    return u32x16_blend<imm>(v, _mm512_sub_epi32(mod2x16(), v));
  }

  static u32x16 reduce(u64x8 x0246, u64x8 x1357) {
    // (x + u64(u32(x) * IR) * MOD) >> 32;
    auto y0246 = u32x16_mul0246(u32x16_mul0246(x0246, irx16()), modx16());
    auto y1357 = u32x16_mul0246(u32x16_mul0246(x1357, irx16()), modx16());
    auto z0246 = _mm512_add_epi64(x0246, y0246);
    z0246 = i32x16_swap_lohi(z0246);
    auto z1357 = _mm512_add_epi64(x1357, y1357);
    return u32x16_blend<0xaaaa>(z0246, z1357);
  }

  static u32x16 mul_reduce(u32x16 a, u32x16 b) {
    // return reduce(u64(a) * b);
    u64x8 x0246 = u32x16_mul0246(a, b);
    u64x8 x1357 = u32x16_mul1357(a, b);
    return reduce(x0246, x1357);
  }

  static u32x16 trans(u32x16 v) {
    return mul_reduce(v, r2x16());
  }

  static u32x16 get(u32x16 v) {
    u32x16 v1 = mul_reduce(v, _mm512_set1_epi32(1));
    u32x16 v2 = _mm512_sub_epi32(v1, modx16());
    return _mm512_min_epu32(v1, v2);
  }

  M32x16 &operator*=(const M32x16 &rhs) {
    v = mul_reduce(v, rhs.v);
    return *this;
  }

  static u32x16 div2(u32x16 v) {
    // (v & 1) ? v + mod : v
    u32x16 one = _mm512_set1_epi32(1);
    u32x16 and1 = _mm512_and_si512(one, v);
    and1 = _mm512_sub_epi32(and1, one);
    u32x16 smodx16 = _mm512_andnot_si512(and1, modx16());
    u32x16 ans = _mm512_add_epi32(v, smodx16);
    return _mm512_srli_epi32(ans, 1);
  }

  static u32x16 irx16() {
    return _mm512_set1_epi32(ModT::IR);
  }
  static u32x16 r2x16() {
    return _mm512_set1_epi32(ModT::R2);
  }
  static u32x16 mod2x16() {
    return _mm512_set1_epi32(ModT::MOD2);
  }
  static u32x16 modx16() {
    return _mm512_set1_epi32(ModT::MOD);
  }

  static void dot(ModT *f, const ModT *g, u32 n) {
    auto *fx = reinterpret_cast<M32x16 *>(f);
    auto *gx = reinterpret_cast<const M32x16 *>(g);
    for (u32 i = 0; i != n / 8; ++i)
      fx[i] *= gx[i];
  }

  static void dot1(ModT *f, u32 n, const ModT g0) {
    auto *fx = reinterpret_cast<M32x16 *>(f);
    auto gx = from(g0);
    for (u32 i = 0; i != n / 8; ++i)
      f[i] *= gx;
  }
};

ALGO_END_NAMESPACE

#endif
