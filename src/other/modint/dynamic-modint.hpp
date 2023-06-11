#ifndef ALGO_MODINT_DYNAMIC_MODINT
#define ALGO_MODINT_DYNAMIC_MODINT

#include "../../base.hpp"
#include "../../math/cipolla.hpp"
#include <iostream>

#ifdef DEBUG
#include <cassert>
#endif

// 封装 Modint，功能由 Space 提供

template <class Space_>
struct DynamicModint {
  using Space = Space_;
  using ValueT = typename Space::ValueT;
  using TransT = typename Space::TransT;
  using isStatic = std::false_type;
  using rawU32 = typename Space::rawU32;
  using isMontgomery = typename Space::isMontgomery;

  TransT v;

  DynamicModint() = default;
  DynamicModint(ValueT v_) : v(Space::trans(v_)) {}

  using Self = DynamicModint;

  static bool set_mod(u32 mod) {
    return Space::set_mod(mod);
  }

  explicit operator ValueT() const {
    return val();
  }

  static Self safe(i64 v) {
    return Space::safe(v);
  }

  static Self raw(u32 v) {
    Self r;
    r.v = v;
    return r;
  }

  ValueT val() const {
    return Space::val(v);
  }

  TransT raw() const {
    return v;
  }

  static ValueT mod() {
    return Space::mod();
  }

  Self &operator+=(const Self &rhs) {
    v = Space::add(v, rhs.v);
    return *this;
  }

  Self &operator-=(const Self &rhs) {
    v = Space::sub(v, rhs.v);
    return *this;
  }

  Self &operator*=(const Self &rhs) {
    v = Space::mul(v, rhs.v);
    return *this;
  }

  friend inline Self operator+(const Self &lhs, const Self &rhs) {
    return Self(lhs) += rhs;
  }

  friend inline Self operator-(const Self &lhs, const Self &rhs) {
    return Self(lhs) -= rhs;
  }

  friend inline Self operator*(const Self &lhs, const Self &rhs) {
    return Self(lhs) *= rhs;
  }

  Self pow(u64 n) const {
    Self r(1), a(*this);
    for (; n > 0; n /= 2) {
      if (n % 2 == 1)
        r *= a;
      a *= a;
    }
    return r;
  }

  Self inv() const {
    return pow(Space::mod() - 2);
  }

  Self &operator/=(const Self &rhs) {
    return *this *= rhs.inv();
  }

  friend inline Self operator/(const Self &lhs, const Self &rhs) {
    return Self(lhs) /= rhs;
  }

  Self operator-() const {
    return Self() -= *this;
  }

  std::optional<Self> sqrt() const {
    return cipola(*this);
  }

  Self shift2() const {
    return Space::shift2(v);
  }

  static Self muladd(const Self &a, const Self &b, const Self &c) { // a * b + c
    return raw(Space::muladd(a.raw(), b.raw(), c.raw()));
  }

  static Self mulsub(const Self &a, const Self &b, const Self &c) { // a * b - c
    return raw(Space::mulsub(a.raw(), b.raw(), c.raw()));
  }

  static Self addmul(const Self &a, const Self &b, const Self &c) { // (a + b) * c
    return raw(Space::addmul(a.raw(), b.raw(), c.raw()));
  }

  static Self submul(const Self &a, const Self &b, const Self &c) { // (a - b) * c
    return raw(Space::submul(a.raw(), b.raw(), c.raw()));
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

#endif // ALGO_MODINT_STATIC_MODINT
