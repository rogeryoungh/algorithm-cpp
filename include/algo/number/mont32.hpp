#pragma once

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

using m32 = u32;

struct Mont32 {
  using U = u32;
  using UU = u64;
  const u32 MOD, MOD2, R, IR, R2, ONE;
	
  explicit constexpr Mont32(u32 mod)
      : MOD(mod), MOD2(mod * 2), R(getR(mod)), IR(-getNR(mod)), R2(u64(R) * R % MOD),
        ONE(trans(1)) {}
  constexpr static u32 getR(u32 mod) {
    return (u64(1) << 32) % mod;
  }
  constexpr static u32 getNR(u32 mod) {
    u32 x = 1;
    for (i32 i = 0; i < 5; ++i)
      x *= 2 - x * mod;
    return x;
  }
  inline constexpr m32 trans(u32 x) const {
    // return (u64(x) << 32) % MOD;
    return reduce(u64(x) * R2);
  }
  inline constexpr m32 reduce(u64 x) const {
    return (x + u64(u32(x) * IR) * MOD) >> 32;
  }
  inline constexpr u32 add(u32 a, u32 b) const {
    u32 v1 = a + b, v2 = v1 - MOD2;
    return i32(v2) < 0 ? v1 : v2;
  }
  inline constexpr u32 sub(u32 a, u32 b) const {
    u32 v1 = a - b, v2 = v1 + MOD2;
    return i32(v2) >= 0 ? v2 : v1;
  }
  inline constexpr m32 mul(m32 a, m32 b) const {
    return reduce(u64(a) * b);
  }
  inline constexpr u32 qpow(u32 a, u64 n, u32 r) const {
    for (; n > 0; n /= 2) {
      if (n % 2 == 1)
        r = mul(r, a);
      a = mul(a, a);
    }
    return r;
  }
  inline constexpr u32 qpow(u32 a, u64 n) const {
    return qpow(a, n, ONE);
  }
  inline constexpr m32 inv(m32 x) const {
    return qpow(x, MOD - 2);
  }
  inline constexpr m32 div(m32 a, m32 b) const {
    return reduce(qpow(b, MOD - 2, a));
  }
  inline constexpr u32 get(m32 x) const {
    u32 v1 = reduce(x), v2 = v1 - MOD;
    return i32(v2) < 0 ? v1 : v2;
  }
  inline constexpr u32 div2(u32 x) const {
    return (x % 2 == 1 ? x + MOD : x) >> 1;
  }
  inline constexpr bool cmp(m32 a, m32 b) const {
    return get(a) == get(b);
  }
  inline constexpr bool ncmp(m32 a, m32 b) const {
		return !cmp(a, b);
  }
  inline constexpr m32 neg(m32 x) const {
    return x != 0 ? MOD2 - x : 0;
  }
};

ALGO_END_NAMESPACE
