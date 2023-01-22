#ifndef ALGO_MATH_BAILLIE_PSW
#define ALGO_MATH_BAILLIE_PSW

#include "../../base.hpp"
#include "miller-rabin.hpp"
#include "lucas.hpp"

bool baillie_psw_test(u64 n) {
  if (n <= 6)
    return n == 2 || n == 3 || n == 5;
  if (n % 6 != 1 && n % 6 != 5)
    return false;
  if (!miller_rabbin_base(n, 2))
    return false;
  u64 d = 5;
  while (jacobi(d, n) != -1)
    d += 4;
  return lucas_test_pd(n, d, d);
}

#endif // ALGO_MATH_BAILLIE_PSW
