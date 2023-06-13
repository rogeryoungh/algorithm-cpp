#ifndef ALGO_MODINT_MONTGOMERY_SPACE_X8
#define ALGO_MODINT_MONTGOMERY_SPACE_X8

#include "../../base.hpp"
#include "../avx2.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <type_traits>

namespace simd {

// 仅在 Montgomery 空间里

template <montgomery_modint_concept ModT_>
struct M32x8 {
  using ModT = ModT_;
  u32x8 v;

  M32x8() = default;

  M32x8(const u32x8 a) : v(a) {}

  template <class U32>
  M32x8(const std::array<U32, 8> &a) {
    static_assert(sizeof(U32) == 4);
    v = *reinterpret_cast<const i256 *>(a.data());
  }

  static M32x8 from(ModT v) {
    return i32x8_set1(v.raw());
  }

  explicit operator u32x8() {
    return v;
  }

  static u32x8 reduce_m(u32x8 v) { // 4, 0.33
    // v = sign(v) ? v + mod : v;
    u32x8 sign = i32x8_sign(v);
    u32x8 smodx8 = _mm256_and_si256(sign, modx8());
    return _mm256_add_epi32(v, smodx8);
  }

  static u32x8 reduce_2m(u32x8 v) { // 4, 0.33
    // v = sign(v) ? v + mod2 : v;
    u32x8 sign = i32x8_sign(v);
    u32x8 smod2x8 = _mm256_and_si256(sign, mod2x8());
    return _mm256_add_epi32(v, smod2x8);
  }

  M32x8 &operator+=(const M32x8 &rhs) { // 6, 0.33
    v = _mm256_add_epi32(v, rhs.v);
    v = _mm256_sub_epi32(v, mod2x8());
    v = reduce_2m(v);
    return *this;
  }

  M32x8 &operator-=(const M32x8 &rhs) { // 5, 0.33
    v = _mm256_sub_epi32(v, rhs.v);
    v = reduce_2m(v);
    return *this;
  }

  M32x8 operator-() const {
    u32x8 sub = _mm256_sub_epi32(mod2x8(), v); // 1, 0.33
    u32x8 sign = i32x8_sign(v);                // 2, 0.33
    return _mm256_andnot_si256(sign, sub);
    // return neg<0b11111111>();
  }

  static u32x8 reduce(const u64x4 &x0246, const u64x4 &x1357) { // 24, 0.33
    // (x + u64(u32(x) * IR) * MOD) >> 32;
    auto y0246 = u32x8_mul0246(u32x8_mul0246(x0246, irx8()), modx8());
    auto y1357 = u32x8_mul0246(u32x8_mul0246(x1357, irx8()), modx8());
    auto z0246 = _mm256_add_epi64(x0246, y0246);
    z0246 = i32x8_swap_lohi(z0246);
    auto z1357 = _mm256_add_epi64(x1357, y1357);
    // z1357 = i32x8::shuffle<0b11110101>(z1357);
    return _mm256_blend_epi32(z0246, z1357, 0b10101010);
  }

  static u32x8 mul_reduce(const u32x8 &a, const u32x8 &b) { // 36, 0.33
    // return reduce(u64(a) * b);
    u64x4 x0246 = u32x8_mul0246(a, b);
    u64x4 x1357 = u32x8_mul1357(a, b);
    return reduce(x0246, x1357);
  }

  // static I32x8 mul_reduce(const I32x8 &a, const U64x4 &b) {
  //   U32x8 x0 = u32x8::mul_lo(irx8(), u32x8::mul_lo(a, b));
  //   auto [x0246, x1357] = u32x8::mul_0246_1357(a, b);
  //   auto x1 = u32x8::mul(x0, modx8());
  //   auto x2 = u32x8::mul(i32x8::shuffle<0b11110101>(x0), modx8());
  //   auto z0246 = u64x4::shift_r<32>(i64x4::add(x0246, x1));
  //   auto z1357 = i64x4::add(x1357, x2);
  //   return i32x8::blend<0b10101010>(z0246, z1357);
  // }

  M32x8 &operator*=(const M32x8 &rhs) {
    v = mul_reduce(v, rhs.v);
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
  static M32x8 addmul(const M32x8 &a, const M32x8 &b, const M32x8 &c) {
    // (a + b) * c
    auto v = _mm256_add_epi32(a.v, b.v);
    return mul_reduce(v, c.v);
    // return (a + b) * c;
  }

  static M32x8 submul(const M32x8 &a, const M32x8 &b, const M32x8 &c) {
    // (a - b) * c
    auto v = _mm256_sub_epi32(a.v, b.v);
    v = _mm256_add_epi32(v, mod2x8());
    return mul_reduce(v, c.v);
    // return (a - b) * c;
  }

  u32x8 raw() const {
    return v;
  }

  template <i32 imm>
  M32x8 neg() const {
    auto sub = _mm256_sub_epi32(mod2x8(), v);
    return _mm256_blend_epi32(v, sub, imm);
  }

  auto to_array() const {
    return i256_toarray<u32>(v);
  }

  template <i32 imm>
  M32x8 shuffle() const {
    return _mm256_shuffle_epi32(v, imm);
  }

  template <i32 imm>
  M32x8 shufflex4() const {
    return _mm256_permute2x128_si256(v, v, imm);
  }

  static u32x8 irx8() {
    return i32x8_set1(ModT::Space::ir());
    // static u32x8 IRx8 = i32x8_set1(ModT::Space::ir());
    // return IRx8;
  }
  static u32x8 mod2x8() {
    return i32x8_set1(ModT::Space::mod2());
    // static u32x8 MOD2x8 = i32x8_set1(ModT::Space::mod2());
    // return MOD2x8;
  }
  static u32x8 modx8() {
    return i32x8_set1(ModT::Space::mod());
    // static u32x8 MODx8 = i32x8_set1(ModT::Space::mod2());
    // return MODx8;
  }
};

} // namespace simd

#endif // ALGO_MODINT_MONTGOMERY_SPACE_X8
