#ifndef ALGO_MATH_GCD_ORIGINAL_REC
#define ALGO_MATH_GCD_ORIGINAL_REC

#include "../../base.hpp"

#include <cassert>
#include <utility>
#include <tuple>

constexpr u64 gcd(u64 x, u64 y) {
  if (x == 0)
    return y;
  else
    return gcd(y % x, x);
}

constexpr u64 lcm(u64 x, u64 y) {
  return x / gcd(x, y) * y;
}

constexpr std::tuple<i64, i64, i64> exgcd(i64 a, i64 b) {
  if (b == 0) {
    if (a > 0)
      return {1, 0, a};
    else
      return {-1, 0, -a};
  } else {
    auto [y, x, g] = exgcd(b, a % b);
    return {x, y - a / b * x, g};
  }
}

#endif // ALGO_MATH_GCD_ORIGINAL_REC
