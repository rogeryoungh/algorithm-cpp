#ifndef ALGO_MATH_GCD_ORIGINAL
#define ALGO_MATH_GCD_ORIGINAL

#include "../../base.hpp"

#include <bit>
#include <cassert>
#include <utility>
#include <tuple>

constexpr u64 gcd(u64 x, u64 y) {
  if (x > y)
    std::swap(x, y);
  while (x != 0)
    std::swap(x, y %= x);
  return y;
}

constexpr u64 lcm(u64 x, u64 y) {
  return x / gcd(x, y) * y;
}

std::tuple<i64, i64, i64> exgcd(i64 a, i64 b);

#endif // ALGO_MATH_GCD_ORIGINAL
