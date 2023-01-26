#ifndef ALGO_MODINT_BASIC_MODINT
#define ALGO_MODINT_BASIC_MODINT

#include "../../base.hpp"
#include "../../math/cipolla.hpp"
#include <iostream>

#ifdef DEBUG
#include <cassert>
#endif

// 朴素的 Modint30，采取了一些常数优化
// 不保证 inv 等运算在 MOD 下合法

template <u32 MOD>
class BasicModint {
protected:
  u32 v;
  enum : u32 { r2 = u64(1) << 32 / MOD };
  static_assert(0 < MOD && MOD < u32(1) << 31, "mod must in [1, 2^31)");
  constexpr inline BasicModint(u32 v_, u32) : v(v_) {}

public:
  using value_type = u32;
  using is_static = std::true_type;
  using raw_u32 = std::true_type;

  constexpr inline static u32 reduce_2m(u32 n) {
#ifdef DEBUG
    assert(n < MOD * 2);
#endif
    return n - (n >= MOD ? MOD : 0);
  }

  constexpr inline static u32 reduce_neg(u32 n) {
    return n + (n >> 31 ? MOD : 0);
  }

  using Self = BasicModint;
  constexpr BasicModint() : Self(0, u32()) {}
  constexpr BasicModint(i64 v_) : BasicModint(reduce_neg(v_ % MOD), u32()) {}

  explicit operator u32() const {
    return val();
  }

  constexpr static Self from_raw(u32 v) {
    return Self(v, u32());
  }

  constexpr u32 val() const {
    return v;
  }

  constexpr static u32 get_mod() {
    return MOD;
  }

  constexpr Self &operator+=(const Self &rhs) {
    v = reduce_neg(v + rhs.v - MOD);
    return *this;
  }

  constexpr Self &operator-=(const Self &rhs) {
    v = reduce_neg(v - rhs.v);
    return *this;
  }

  constexpr Self &operator*=(const Self &rhs) {
    v = u64(v) * rhs.v % MOD;
    return *this;
  }

  constexpr Self operator+(const Self &b) const {
    return Self(*this) += b;
  }

  constexpr Self operator-(const Self &b) const {
    return Self(*this) -= b;
  }

  constexpr Self operator*(const Self &b) const {
    return Self(*this) *= b;
  }

  constexpr Self pow(u64 n) const {
    Self r(1), a(*this);
    for (; n > 0; n /= 2) {
      if (n % 2 == 1)
        r *= a;
      a *= a;
    }
    return r;
  }

  constexpr Self inv() const {
    return pow(MOD - 2);
  }

  constexpr Self &operator/=(const Self &rhs) {
    return *this *= rhs.inv();
  }

  constexpr Self operator/(const Self &b) const {
    return Self(*this) /= b;
  }

  constexpr Self operator-() const {
    return Self{v == 0 ? 0 : MOD - v};
  }

  constexpr std::optional<Self> sqrt() const {
    return cipola(val(), get_mod());
  }

  constexpr Self shift2() const {
    return (v + (v % 2 == 0 ? 0 : get_mod())) / 2;
  }

  friend inline std::istream &operator>>(std::istream &is, Self &m) {
    i64 x;
    is >> x;
    m = Self(x);
    return is;
  }

  friend inline std::ostream &operator<<(std::ostream &os, const Self &m) {
    return os << m.val();
  }

  friend inline bool operator==(const Self &lhs, const Self &rhs) {
    return lhs.val() == rhs.val();
  }

  friend inline bool operator!=(const Self &lhs, const Self &rhs) {
    return !(lhs == rhs);
  }
};

#include "modint-io.hpp"

#endif // ALGO_MODINT_BASIC_MODINT
