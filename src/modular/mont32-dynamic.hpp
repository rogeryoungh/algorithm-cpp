#ifndef ALGO_H_MODULAR_MONT32_DYNAMIC
#define ALGO_H_MODULAR_MONT32_DYNAMIC

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <i32>
struct M32D {
  inline static u32 MOD, MOD2, R, IR, R2;
  u32 v;

  static bool setMod(u32 mod) {
    if (mod <= 1 || mod >= (1u << 30))
      return false;
    if (mod % 2 == 0)
      return false;
    MOD = mod, MOD2 = MOD * 2;
    R = (u64(1) << 32) % MOD;
    R2 = (u64(R) * R) % MOD;
    IR = -getNR();
    return true;
  }

  M32D(u32 val = 0) : v(trans(val)) {}

  static M32D from_raw(u32 x) {
    return reinterpret_cast<M32D>(x);
  }

  u32 raw() const {
    return v;
  }

  static u32 getNR() {
    u32 x = 1;
    for (i32 i = 0; i < 5; ++i)
      x *= 2 - x * MOD;
    return x;
  }

  static u32 trans(u32 x) {
    return reduce(u64(x) * R2);
  }

  static u32 reduce(u64 x) {
    return (x + u64(u32(x) * IR) * MOD) >> 32;
  }

  static u32 reduce_m(u32 n) {
    return n >> 31 ? n + MOD : n;
  }

  static u32 reduce_2m(u32 n) {
    return n >> 31 ? n + MOD2 : n;
  }

  u32 get() const {
    return reduce_m(reduce(v) - MOD);
  }

  M32D &operator+=(M32D o) {
    v = reduce_2m(v + o.v - MOD2);
    return *this;
  }

  M32D &operator-=(M32D o) {
    v = reduce_2m(v - o.v);
    return *this;
  }

  M32D &operator*=(M32D o) {
    v = reduce(u64(v) * o.v);
    return *this;
  }

  M32D operator-() const {
    return from_raw(v != 0 ? MOD2 - v : v);
  }

  M32D div2() const {
    return from_raw((v % 2 == 1 ? v + MOD : v) >> 1);
  }

  M32D pow(u64 n, M32D r = {1}) const {
    M32D a = *this;
    for (; n > 0; n /= 2) {
      if (n % 2 == 1)
        r *= a;
      a *= a;
    }
    return r;
  }

  M32D inv() const {
    return pow(MOD - 2);
  }

  M32D &operator/=(M32D o) {
    v = pow(MOD - 2, o);
    return *this;
  }

  bool operator==(M32D o) const {
    return get() == o.get();
  }

  bool operator!=(M32D o) const {
    return get() != o.get();
  }

  friend M32D operator+(const M32D &l, const M32D &r) {
    return M32D(l) += r;
  }
  friend M32D operator-(const M32D &l, const M32D &r) {
    return M32D(l) -= r;
  }
  friend M32D operator*(const M32D &l, const M32D &r) {
    return M32D(l) *= r;
  }
  friend M32D operator/(const M32D &l, const M32D &r) {
    return M32D(l) /= r;
  }
};

ALGO_END_NAMESPACE

#endif
