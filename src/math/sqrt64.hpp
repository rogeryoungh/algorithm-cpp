#ifndef ALGO_MATH_SQRT64
#define ALGO_MATH_SQRT64

#include "../base.hpp"

#include <utility>

u32 sqrt64(u64 n) {
  u64 x = n, x2 = n;
  do {
    x2 = (x + n / x) / 2;
    std::swap(x, x2);
  } while (x < x2);
  return x;
}

#endif // ALGO_MATH_SQRT64
