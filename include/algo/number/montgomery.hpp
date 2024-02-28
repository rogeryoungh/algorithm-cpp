#pragma once

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <class U>
struct Mont {
  using S = std::make_signed_t<U>;
  using UU = std::conditional_t<std::is_same_v<U, u32>, u64, u128>;
  const U MOD, MOD2, R, IR, R2, ONE;

  explicit constexpr Mont(U mod)
      : MOD(mod), MOD2(mod * 2), R(getR(mod)), IR(-getNR(mod)), R2(UU(R) * R % MOD), ONE(trans(1)) {
  }
  constexpr static U getR(U mod) {
    return (UU(1) << (sizeof(U) * 8)) % mod;
  }
  constexpr static U getNR(U mod) {
    U x = 1;
    for (u32 i = 0; i != 6; ++i)
      x *= 2 - x * mod;
    return x;
  }
  constexpr U trans(U x) const {
    // return (u64(x) << 32) % MOD;
    return reduce(UU(x) * R2);
  }
  constexpr U reduce(UU x) const {
    return (x + UU(U(x) * IR) * MOD) >> (sizeof(U) * 8);
  }
  constexpr U norm(U v) const {
    U v2 = v - MOD;
    return S(v2) < 0 ? v : v2;
  }
  constexpr U add(U a, U b) const {
    U v1 = a + b, v2 = v1 - MOD2;
    return S(v2) < 0 ? v1 : v2;
  }
  constexpr U sub(U a, U b) const {
    U v1 = a - b, v2 = v1 + MOD2;
    return S(v1) >= 0 ? v1 : v2;
  }
  constexpr U mul(U a, U b) const {
    return reduce(UU(a) * b);
  }
  constexpr U qpow(U a, u64 n, U r) const {
    for (; n > 0; n /= 2) {
      if (n % 2 == 1)
        r = mul(r, a);
      a = mul(a, a);
    }
    return r;
  }
  constexpr U qpow(U a, u64 n) const {
    return qpow(a, n, ONE);
  }
  constexpr U inv(U x) const {
    return qpow(x, MOD - 2);
  }
  constexpr U div(U a, U b) const {
    return reduce(qpow(b, MOD - 2, a));
  }
  constexpr U get(U x) const {
    return norm(reduce(x));
  }
  constexpr U div2(U x) const {
    return (x % 2 == 1 ? x + MOD : x) >> 1;
  }
  constexpr bool cmp(U a, U b) const {
    return get(a) == get(b);
  }
  constexpr bool ncmp(U a, U b) const {
    return !cmp(a, b);
  }
  constexpr U neg(U x) const {
    return x != 0 ? MOD2 - x : 0;
  }
};

using Mont32 = Mont<u32>;
using Mont64 = Mont<u64>;

template <class ModT>
using ModU = ModT::U;

template <class ModT>
using ModUU = ModT::UU;

ALGO_END_NAMESPACE
