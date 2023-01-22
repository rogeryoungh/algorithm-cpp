#ifndef ALGO_MATH_QPOW_U32
#define ALGO_MATH_QPOW_U32

#include "../../base.hpp"

constexpr u32 qpow(u32 a, u64 b, u32 m) {
  u32 r = 1;
  for (; b > 0; b /= 2) {
    if (b % 2 == 1)
      r = u64(a) * r % m;
    a = u64(a) * a % m;
  }
  return r;
}

#endif // ALGO_MATH_QPOW_U32
