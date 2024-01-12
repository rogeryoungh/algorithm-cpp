#ifndef ALGO_H_MODULAR_MONT32_CONST
#define ALGO_H_MODULAR_MONT32_CONST

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <u32 _M>
struct M32C {

  inline static constexpr u32 MOD = _M, MOD2 = MOD * 2, R = (u64(1) << 32) % MOD;
  static constexpr u32 get_nr() {
    u32 x = 1;
    for (i32 i = 0; i < 5; ++i)
      x *= 2 - x * MOD;
    return x;
  }
  inline static constexpr u32 R2 = u64(R) * R % MOD, IR = -get_nr();

  u32 v;

  constexpr M32C(u32 val = 0) : v(trans(val)) {}

  static constexpr M32C from_raw(u32 x) {
    return std::bit_cast<M32C>(x);
  }

  constexpr u32 raw() const {
    return v;
  }

  static constexpr u32 trans(u32 x) {
    return reduce(u64(x) * R2);
  }

  static constexpr u32 get(u32 x) {
    return reduce_m(reduce(x) - MOD);
  }

  static constexpr u32 reduce(u64 x) {
    return (x + u64(u32(x) * IR) * MOD) >> 32;
  }

  static constexpr u32 reduce_m(u32 n) {
    return n >> 31 ? n + MOD : n;
  }

  static constexpr u32 reduce_2m(u32 n) {
    return n >> 31 ? n + MOD2 : n;
  }

  constexpr u32 get() const {
    return reduce_m(reduce(v) - MOD);
  }

  constexpr M32C &operator+=(M32C o) {
    u32 v1 = v + o.v, v2 = v1 - MOD2;
    v = i32(v2) < 0 ? v1 : v2;
    return *this;
  }

  constexpr M32C &operator-=(M32C o) {
    u32 v1 = v - o.v, v2 = v1 + MOD2;
    v = i32(v1) < 0 ? v2 : v1;
    return *this;
  }

  constexpr M32C &operator*=(M32C o) {
    v = reduce(u64(v) * o.v);
    return *this;
  }

  constexpr M32C operator-() const {
    return from_raw(v != 0 ? MOD2 - v : v);
  }

  constexpr M32C div2() const {
    return from_raw((v % 2 == 1 ? v + MOD : v) >> 1);
  }

  constexpr M32C pow(u64 n, M32C r = {1}) const {
    M32C a = *this;
    for (; n > 0; n /= 2) {
      if (n % 2 == 1)
        r *= a;
      a *= a;
    }
    return r;
  }

  constexpr M32C inv() const {
    return pow(MOD - 2);
  }

  constexpr M32C &operator/=(M32C o) {
    v = pow(MOD - 2, o);
    return *this;
  }

  constexpr bool operator==(M32C o) const {
    return get() == o.get();
  }

  constexpr bool operator!=(M32C o) const {
    return get() != o.get();
  }

  friend constexpr M32C operator+(const M32C &l, const M32C &r) {
    return M32C(l) += r;
  }
  friend constexpr M32C operator-(const M32C &l, const M32C &r) {
    return M32C(l) -= r;
  }
  friend constexpr M32C operator*(const M32C &l, const M32C &r) {
    return M32C(l) *= r;
  }
  friend constexpr M32C operator/(const M32C &l, const M32C &r) {
    return M32C(l) /= r;
  }
};

ALGO_END_NAMESPACE

#endif
