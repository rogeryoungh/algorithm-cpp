#ifndef ALGO_MATH_MILLER_RABIN
#define ALGO_MATH_MILLER_RABIN

#include "../../base.hpp"
#include "../qpow/u64.hpp"

#include <bit>
#include <array>

bool miller_rabbin_base(u64 n, u64 a) {
  u32 s = std::countr_zero(n - 1);
  u64 d = n >> s;
  u64 ad = qpow64(a, d, n);
  if (ad == 1 || ad == n - 1 || ad == 0)
    return true;
  for (u32 i = 1; i <= s - 1; ++i) {
    ad = u128(ad) * ad % n;
    if (ad == n - 1)
      return true;
    if (ad == 1)
      break;
  }
  return false;
}

bool miller_rabbin_test(u64 n) {
  if (n <= 6)
    return n == 2 || n == 3 || n == 5;
  if (n % 6 != 1 && n % 6 != 5)
    return false;
  for (u32 a : {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {
    if (!miller_rabbin_base(n, a))
      return false;
  }
  return true;
}

#endif // ALGO_MATH_MILLER_RABIN
