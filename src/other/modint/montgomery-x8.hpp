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
  I32x8 v;

  M32x8() = default;

  M32x8(const I32x8 &a) : v(a) {}

  template <class U32>
  M32x8(const std::array<U32, 8> &a) {
    static_assert(sizeof(U32) == 4);
    v = i256::load((const I256 *)a.data());
  }

  static M32x8 from(i32 v) {
    return i32x8::from(v);
  }

  static M32x8 from(ModT v) {
    return from(v.raw());
  }

  explicit operator I32x8() {
    return v;
  }

  static I32x8 get_irx8() {
    return i32x8::from(ModT::Space::ir());
  }

  static I32x8 get_mod2x8() {
    return i32x8::from(ModT::Space::mod2());
  }

  static I32x8 get_modx8() {
    return i32x8::from(ModT::Space::mod());
  }

  static I32x8 reduce_m(I32x8 v) {
    I32x8 sign = i32x8::sign(v);
    v = i32x8::add(v, i256::bit_and(sign, get_modx8()));
    return v;
  }

  static I32x8 reduce_2m(I32x8 v) {
    I32x8 sign = i32x8::sign(v);
    v = i32x8::add(v, i256::bit_and(sign, get_mod2x8()));
    return v;
  }

  M32x8 &operator+=(const M32x8 &rhs) {
    v = i32x8::add(v, rhs.v);
    v = i32x8::sub(v, get_mod2x8());
    v = reduce_2m(v);
    return *this;
  }

  M32x8 &operator-=(const M32x8 &rhs) {
    v = i32x8::sub(v, rhs.v);
    v = reduce_2m(v);
    return *this;
  }

  M32x8 operator-() const {
    return neg<0b11111111>();
  }

  static I32x8 reduce(const U64x4 &x0246, const U64x4 &x1357) {
    auto km0246 = u32x8::mul(u32x8::mul(x0246, get_irx8()), get_modx8());
    auto km1357 = u32x8::mul(u32x8::mul(x1357, get_irx8()), get_modx8());
    auto z0246 = i64x4::add(x0246, km0246);
    z0246 = i32x8::shuffle<0b11110101>(z0246);
    auto z1357 = i64x4::add(x1357, km1357);
    // z1357 = i32x8::shuffle<0b11110101>(z1357);
    return i32x8::blend<0b10101010>(z0246, z1357);
  }

  M32x8 &operator*=(const M32x8 &rhs) {
    auto [x0246, x1357] = u32x8::mul_0246_1357(v, rhs.v);
    v = reduce(x0246, x1357);
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
  constexpr static M32x8 addmul(const M32x8 &a, const M32x8 &b, const M32x8 &c) { // (a + b) * c
    auto v = i32x8::add(a.v, b.v);
    auto [x0246, x1357] = u32x8::mul_0246_1357(v, c.v);
    v = reduce(x0246, x1357);
    return v;
  }

  constexpr static M32x8 submul(const M32x8 &a, const M32x8 &b, const M32x8 &c) { // (a - b) * c
    auto v = i32x8::sub(a.v, b.v);
    v = i32x8::add(v, get_mod2x8());
    auto [x0246, x1357] = u32x8::mul_0246_1357(v, c.v);
    v = reduce(x0246, x1357);
    return v;
  }

  U32x8 raw() const {
    return v;
  }

  template <i32 imm>
  M32x8 neg() const {
    auto sub = i32x8::sub(get_mod2x8(), v);
    return i32x8::blend<imm>(v, sub);
  }

  auto to_array() const {
    return i256::to_array<u32>(v);
  }

  template <i32 imm>
  M32x8 shuffle() const {
    return i32x8::shuffle<imm>(v);
  }

  template <i32 imm>
  M32x8 shufflex4() const {
    return i128x2::shuffle<imm>(v);
  }
};

} // namespace simd

#endif // ALGO_MODINT_MONTGOMERY_SPACE_X8
