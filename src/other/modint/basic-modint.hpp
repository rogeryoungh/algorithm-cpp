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
  static_assert(0 < MOD && MOD < u32(1) << 31, "mod must in [1, 2^31)");

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
  constexpr BasicModint(u32 v_ = 0) : v(v_) {}

  explicit operator u32() const {
    return val();
  }

  constexpr static Self safe(i64 v) {
    return Self(reduce_neg(v % MOD));
  }

  template <class U = u32>
  constexpr inline U val() const {
    return U(v);
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

  friend constexpr inline Self operator+(const Self &lhs, const Self &rhs) {
    return Self(lhs) += rhs;
  }

  friend constexpr inline Self operator-(const Self &lhs, const Self &rhs) {
    return Self(lhs) -= rhs;
  }

  friend constexpr inline Self operator*(const Self &lhs, const Self &rhs) {
    return Self(lhs) *= rhs;
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

  friend constexpr inline Self operator/(const Self &lhs, const Self &rhs) {
    return Self(lhs) /= rhs;
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
    m = Self::safe(x);
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
