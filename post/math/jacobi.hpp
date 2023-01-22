#ifndef ALGO_MATH_JACOBI
#define ALGO_MATH_JACOBI

#include "../base.hpp"

#include <cassert>
#include <utility>

// a/n is represented as (a,n)
// https://en.wikipedia.org/wiki/Jacobi_symbol
i32 jacobi(u64 a, u64 n) {
  assert(n > 0 && n % 2 == 1);
  a %= n;
  i32 t = 1;
  while (a != 0) {
    while (a % 2 == 0) {
      a /= 2;
      u64 r = n % 8;
      if (r == 3 || r == 5)
        t = -t;
    }
    std::swap(a, n);
    if (a % 4 == 3 && n % 4 == 3)
      t = -t;
    a %= n;
  }
  return n == 1 ? t : 0;
}

#endif // ALGO_MATH_JACOBI
