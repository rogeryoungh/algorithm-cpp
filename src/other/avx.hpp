#ifndef ALGO_AVX
#define ALGO_AVX

#include "../base.hpp"

// https://judge.yosupo.jp/submission/92714

#include <immintrin.h>
#include <array>

constexpr u32 M998 = 998244353;

using u8x32 = __m256i;
struct i32x8 {
  u8x32 v;
  constexpr i32x8() : v{} {};
  i32x8(u32 x) : v(_mm256_set1_epi32(x)) {}
  constexpr i32x8(const u8x32 &x) : v(x) {}
  template <class U32>
  i32x8(const std::array<U32, 8> &x) : v(_mm256_loadu_si256((u8x32 *)x.data())) {}
  i32x8(const u8x32 *p) : v(_mm256_loadu_si256(p)) {}
  i32x8 &operator+=(const i32x8 rhs) {
    v = _mm256_add_epi32(v, rhs.v);
    return *this;
  }
  i32x8 &operator-=(const i32x8 rhs) {
    v = _mm256_sub_epi32(v, rhs.v);
    return *this;
  }
  i32x8 &operator*=(const i32x8 rhs) {
    v = _mm256_mul_epi32(v, rhs.v);
    return *this;
  }
  i32x8 &operator&=(const i32x8 rhs) {
    v = _mm256_and_si256(v, rhs.v);
    return *this;
  }
  i32x8 operator>(const i32x8 &rhs) const {
    return _mm256_cmpgt_epi32(v, rhs.v);
  }
  friend i32x8 operator+(const i32x8 lhs, const i32x8 rhs) {
    return i32x8(lhs) += rhs;
  }
  friend i32x8 operator-(const i32x8 lhs, const i32x8 rhs) {
    return i32x8(lhs) -= rhs;
  }
  friend i32x8 operator*(const i32x8 lhs, const i32x8 rhs) {
    return i32x8(lhs) *= rhs;
  }
  friend i32x8 operator&(const i32x8 lhs, const i32x8 rhs) {
    return i32x8(lhs) &= rhs;
  }
  i32x8 operator<(const i32x8 &rhs) const {
    return rhs > *this;
  }
  auto to_arr() const {
    alignas(32) std::array<u32, 8> b;
    _mm256_store_si256((__m256i *)b.data(), v);
    return b;
  }
  i32x8 sign() const {
    return *this < i32x8(_mm256_setzero_si256());
  }
  template <u32 imm>
  i32x8 shuffle() const {
    return _mm256_shuffle_epi32(v, imm);
  }
  static std::pair<u8x32, u8x32> mul_0246_1357(const i32x8 lhs, const i32x8 rhs) {
    i32x8 x0246 = lhs * rhs;
    // 76543210 -> 77553311
    i32x8 x1357 = lhs.shuffle<0b11110101>() * rhs.shuffle<0b11110101>();
    return {(u8x32 &)x0246, (u8x32 &)x1357};
  }
  void store(u8x32 *p) {
    _mm256_storeu_si256(p, v);
  }
  template <u32 imm>
  static i32x8 blend(const i32x8 lhs, const i32x8 rhs) {
    return _mm256_blend_epi32(lhs.v, rhs.v, imm);
  }
  template <u32 imm>
  i32x8 permute_64() {
    return _mm256_permute4x64_epi64(v, imm);
  }
  template <u32 imm>
  static i32x8 permute_128(const i32x8 lhs, const i32x8 rhs) {
    return _mm256_permute2x128_si256(lhs.v, rhs.v, imm);
  }
  static i32x8 unpack_hi_64(const i32x8 lhs, const i32x8 rhs) {
    return _mm256_unpackhi_epi64(lhs.v, rhs.v);
  }
  static i32x8 unpack_lo_64(const i32x8 lhs, const i32x8 rhs) {
    return _mm256_unpacklo_epi64(lhs.v, rhs.v);
  }
  static i32x8 unpack_hi_32(const i32x8 lhs, const i32x8 rhs) {
    return _mm256_unpackhi_epi32(lhs.v, rhs.v);
  }
  static i32x8 unpack_lo_32(const i32x8 lhs, const i32x8 rhs) {
    return _mm256_unpacklo_epi32(lhs.v, rhs.v);
  }
};

struct u32x8 : i32x8 {
  using i32x8::i32x8;
  u32x8 &operator*=(const u32x8 &rhs) {
    v = _mm256_mul_epu32(v, rhs.v);
    return *this;
  }
  friend u32x8 operator*(const u32x8 &lhs, const u32x8 &rhs) {
    return u32x8(lhs) *= rhs;
  }
};

struct m32x8 {
  i32x8 v;
  m32x8(const i32x8 &v_) : v(v_) {}
  enum : u32 { P = 998244353, MR = 3296722945 };
  m32x8 &operator+=(const m32x8 &rhs) {
    v += rhs.v;
    v += v.sign() & i32x8(P * 2);
    v -= i32x8(P);
    return *this;
  }
  m32x8 &operator-=(const m32x8 rhs) {
    v -= rhs.v;
    v += v.sign() & i32x8(P * 2);
    v -= i32x8(P);
    return *this;
  }
  m32x8 &operator*=(const m32x8 rhs) {
    i32x8 a = v + i32x8(P);
    i32x8 b = rhs.v + i32x8(P);
    auto [x0246_, x1357_] = i32x8::mul_0246_1357(a, b);
    auto x0246 = (u32x8 &)x0246_;
    auto x1357 = (u32x8 &)x1357_;
    auto km0246 = x0246 * u32x8(MR) * u32x8(P);
    auto km1357 = x1357 * u32x8(MR) * u32x8(P);
    auto z0246 = x0246 - km0246;
    auto z1357 = x1357 - km1357;
    v = i32x8::blend<0b10101010>(z0246.shuffle<0b11110101>(), z1357);
    return *this;
  }
  friend m32x8 operator+(const m32x8 lhs, const m32x8 rhs) {
    return m32x8(lhs) += rhs;
  }
  friend m32x8 operator-(const m32x8 lhs, const m32x8 rhs) {
    return m32x8(lhs) -= rhs;
  }
  friend m32x8 operator*(const m32x8 lhs, const m32x8 rhs) {
    return m32x8(lhs) *= rhs;
  }
};

#endif
