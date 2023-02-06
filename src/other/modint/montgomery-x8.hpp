#ifndef ALGO_MODINT_MONTGOMERY_SPACE_X8
#define ALGO_MODINT_MONTGOMERY_SPACE_X8

#include "../../base.hpp"
#include <type_traits>

#include "../avx2.hpp"

#include "montgomery-space.hpp"

namespace simd {

// 仅在 Montgomery 空间里

template <class ModT>
struct M32x8 {
  I32x8 v;

  M32x8() : v() {}

  M32x8(const I32x8 &a) : v(a) {}

  template <class S>
  M32x8(const M32x8<S> &a) : v(a.v) {}

  template <class U32>
  M32x8(const std::array<U32, 8> &a) {
    static_assert(sizeof(U32) == 4);
    v = i256::load((const I256 *)a.data());
  }

  template <bool aligned = false>
  static M32x8 load(const I256 *p) {
    M32x8 r;
    if constexpr (aligned) {
      r = i256::load(p);
    } else {
      r = i256::loadu(p);
    }
    return r;
  }

  static M32x8 from(i32 v) {
    return i32x8::from(v);
  }

  static M32x8 from(ModT v) {
    return from(v.raw());
  }

  inline static I32x8 Rx8 = i32x8::from(ModT::Space::R);
  inline static I32x8 IRx8 = i32x8::from(ModT::Space::IR);
  inline static I32x8 MOD2x8 = i32x8::from(ModT::Space::MOD2);
  inline static I32x8 MODx8 = i32x8::from(ModT::Space::mod());

  M32x8 &operator+=(const M32x8 &rhs) {
    v = i32x8::add(v, rhs.v);
    v = i32x8::sub(v, MOD2x8);
    I32x8 sign = i32x8::sign(v);
    v = i32x8::add(v, i256::bit_and(sign, MOD2x8));
    return *this;
  }

  M32x8 &operator-=(const M32x8 &rhs) {
    v = i32x8::sub(v, rhs.v);
    I32x8 sign = i32x8::sign(v);
    v = i32x8::add(v, i256::bit_and(sign, MOD2x8));
    return *this;
  }

  static I32x8 reduce(const I64x4 &x0246, const I64x4 &x1357) {
    auto km0246 = u32x8::mul(u32x8::mul(x0246, IRx8), MODx8);
    auto km1357 = u32x8::mul(u32x8::mul(x1357, IRx8), MODx8);
    auto z0246 = i64x4::add(x0246, km0246);
    z0246 = i32x8::shuffle<0b11110101>(z0246);
    auto z1357 = i64x4::add(x1357, km1357);
    z1357 = i32x8::shuffle<0b11110101>(z1357);
    return i32x8::blend<0b10101010>(z0246, z1357);
  }

  M32x8 &operator*=(const M32x8 &rhs) {
    auto [x0246, x1357] = i32x8::mul_0246_1357(v, rhs.v);
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

  I32x8 raw() const {
    return v;
  }

  template <i32 imm>
  M32x8 neg() const {
    auto m2 = i32x8::blend<imm>(i32x8::zero(), MOD2x8);
    return i32x8::abs(i32x8::sub(v, m2));
  }

  template <bool aligned = false>
  void store(I256 *p) {
    if constexpr (aligned) {
      i256::store(p, v);
    } else {
      i256::storeu(p, v);
    }
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
