#pragma once

#include "../number/montgomery.hpp"

#include <bit>
#include <array>

ALGO_BEGIN_NAMESPACE

inline bool basic_primality_test(u64 n) {
  if (n < 64)
    return (0x28208a20a08a28ac >> n) & 1;
  if ((0xdf75d77d >> (n % 30)) & 1)
    return false;
  if (n % 7 == 0 || n % 11 == 0)
    return false;
  return true;
}

template <class U>
inline bool miller_rabin_base(const Mont<U> &_M, TI<U> a, u32 s, TI<U> n) {
  const Mont<U> M = _M;
  U ad = M.qpow(a, n >> s), ad0 = M.get(ad);
  if (ad0 == 1 || ad0 == n - 1 || ad0 == 0)
    return true;
  for (u32 i = 1; i != s; ++i) {
    ad = M.mul(ad, ad), ad0 = M.get(ad);
    if (ad0 == n - 1)
      return true;
    if (ad0 == 1)
      break;
  }
  return false;
}

bool miller_rabin(u64 n) {
  if (!basic_primality_test(n))
    return false;
  u32 s = std::countr_zero(n - 1);
  constexpr std::array bases32 = {2, 7, 61};
  constexpr std::array bases64 = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  if (n < (1 << 30)) {
    auto M = Mont32{u32(n)};
    for (u32 a : bases32) {
      if (!miller_rabin_base(M, a, s, n))
        return false;
    }
  } else {
    auto M = Mont64{n};
    for (u64 a : bases64) {
      if (!miller_rabin_base(M, a, s, n))
        return false;
    }
  }
  return true;
}

ALGO_END_NAMESPACE
