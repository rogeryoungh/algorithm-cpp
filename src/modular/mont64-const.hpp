#ifndef ALGO_H_MODULAR_MONT64_CONST
#define ALGO_H_MODULAR_MONT64_CONST

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <u64 _M>
struct M64C {

  inline static constexpr u64 MOD = _M, MOD2 = MOD * 2, R = (u128(1) << 64) % MOD;
  static constexpr u64 getNR() {
    u64 x = 1;
    for (i64 i = 0; i < 7; ++i)
      x *= 2 - x * MOD;
    return x;
  }
  inline static constexpr u64 R2 = u128(R) * R % MOD, IR = -getNR();

  u64 v;

  M64C(u64 val = 0) : v(trans(val)) {}

  static constexpr M64C from_raw(u64 x) {
    return std::bit_cast<M64C>(x);
  }

  constexpr u64 raw() const {
    return v;
  }

  static constexpr u64 trans(u64 x) {
    return reduce(u128(x) * R2);
  }

  static constexpr u64 get(u64 x) {
    return reduce_m(reduce(x) - MOD);
  }

  static constexpr u64 reduce(u128 x) {
    return (x + u128(u64(x) * IR) * MOD) >> 64;
  }

  static constexpr u64 reduce_m(u64 n) {
    return n >> 63 ? n + MOD : n;
  }

  static constexpr u64 reduce_2m(u64 n) {
    return n >> 63 ? n + MOD2 : n;
  }

  constexpr u64 get() const {
    return reduce_m(reduce(v) - MOD);
  }

  constexpr M64C &operator+=(M64C o) {
    u64 v1 = v + o.v, v2 = v1 - MOD2;
    v = i64(v2) < 0 ? v1 : v2;
    return *this;
  }

  constexpr M64C &operator-=(M64C o) {
    u64 v1 = v - o.v, v2 = v1 + MOD2;
    v = i64(v1) < 0 ? v2 : v1;
    return *this;
  }

  constexpr M64C &operator*=(M64C o) {
    v = reduce(u128(v) * o.v);
    return *this;
  }

  constexpr M64C operator-() const {
    return from_raw(v != 0 ? MOD2 - v : v);
  }

  constexpr M64C div2() const {
    return from_raw((v % 2 == 1 ? v + MOD : v) >> 1);
  }

  constexpr M64C pow(u64 n, M64C r = {1}) const {
    M64C a = *this;
    for (; n > 0; n /= 2) {
      if (n % 2 == 1)
        r *= a;
      a *= a;
    }
    return r;
  }

  constexpr M64C inv() const {
    return pow(MOD - 2);
  }

  constexpr M64C &operator/=(M64C o) {
    v = pow(MOD - 2, o);
    return *this;
  }

  constexpr bool operator==(M64C o) const {
    return get() == o.get();
  }

  constexpr bool operator!=(M64C o) const {
    return get() != o.get();
  }

  friend M64C operator+(const M64C &l, const M64C &r) {
    return M64C(l) += r;
  }
  friend M64C operator-(const M64C &l, const M64C &r) {
    return M64C(l) -= r;
  }
  friend M64C operator*(const M64C &l, const M64C &r) {
    return M64C(l) *= r;
  }
  friend M64C operator/(const M64C &l, const M64C &r) {
    return M64C(l) /= r;
  }
};

ALGO_END_NAMESPACE

#endif
