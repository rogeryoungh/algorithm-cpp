#ifndef ALGO_MATH_LUCAS_TEST
#define ALGO_MATH_LUCAS_TEST

#include "../../base.hpp"
#include "../qpow/u64.hpp"
#include "../jacobi.hpp"
#include "../square-test.hpp"

#include <functional>

bool lucas_test_pd(u64 N, u64 P, u64 D) { // jacobi(D, N) == -1
  using pii = std::pair<u64, u64>;
  auto shift2 = [N](u64 n) {
    return n % 2 == 0 ? n / 2 : (N + n) / 2;
  };
  std::function<pii(u64)> calc = [&](u64 n) {
    if (n == 1)
      return pii{P, 1};
    auto [v, u] = calc(n / 2);
    u64 u2 = u128(u) * v % N;
    u64 v2 = (u128(v) * v + u128(u) * u % N * D) % N;
    v2 = shift2(v2);
    if (n % 2 == 1) {
      u64 u21 = (u128(P) * u2 + v2) % N;
      u64 v21 = (u128(D) * u2 + u128(P) * v2) % N;
      u2 = shift2(u21), v2 = shift2(v21);
    }
    return pii{v2, u2};
  };
  auto [v, u] = calc(N + 1);
  return u == 0;
}

bool lucas_test(u64 n) {
  if (n <= 6)
    return n == 2 || n == 3 || n == 5;
  if (n % 6 != 1 && n % 6 != 5)
    return false;
  if (square_test(n))
    return false;
  for (u64 d = 5, c = 0; d < n && c < 4; d += 4) {
    i32 j = jacobi(d, n);
    if (j == 0)
      return false;
    if (j == 1)
      continue;
    if (!lucas_test_pd(n, d, d))
      return false;
    c++;
  }
  return true;
}

#endif // ALGO_MATH_LUCAS_TEST
