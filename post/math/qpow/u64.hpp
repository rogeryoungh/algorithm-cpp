#ifndef ALGO_MATH_QPOW_U64
#define ALGO_MATH_QPOW_U64

#include "../../base.hpp"

constexpr u64 qpow64(u64 a, u64 b, u64 m) {
  u64 r = 1 % m;
  for (; b > 0; b /= 2) {
    if (b % 2 == 1)
      r = u128(a) * r % m;
    a = u128(a) * a % m;
  }
  return r;
}

#endif // ALGO_MATH_QPOW_U64
