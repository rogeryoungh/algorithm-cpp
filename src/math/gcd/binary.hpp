#ifndef ALGO_MATH_GCD_BINARY
#define ALGO_MATH_GCD_BINARY

#include "../../base.hpp"

#include <algorithm>
#include <bit>
#include <utility>
#include <tuple>

constexpr u64 gcd(u64 x, u64 y) {
  if (x == 0) {
    return y;
  } else if (y == 0) {
    return x;
  } else {
    u32 kx = std::countr_zero(x);
    u32 ky = std::countr_zero(y);
    x >>= kx;
    while (y != 0) {
      y >>= std::countr_zero(y);
      if (x > y)
        std::swap(x, y);
      y -= x;
    }
    return x << std::min(kx, ky);
  }
}

constexpr u64 lcm(u64 x, u64 y) {
  return x / gcd(x, y) * y;
}

std::tuple<i64, i64, i64> exgcd(i64 a, i64 b);

#endif // ALGO_MATH_GCD_BINARY
