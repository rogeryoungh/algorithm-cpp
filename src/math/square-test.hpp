#ifndef ALGO_MATH_SQUARE_TEST
#define ALGO_MATH_SQUARE_TEST

#include "../base.hpp"
#include "sqrt64.hpp"

#include <array>

namespace detail {

template <u32 u>
constexpr auto gen_square_table() {
  std::array<bool, u> arr{};
  for (u32 i = 0; i < u; i++)
    arr[i * i % u] = 1;
  return arr;
}

} // namespace detail

bool square_test(u64 n) {
  static auto q11 = detail::gen_square_table<11>();
  static auto q63 = detail::gen_square_table<63>();
  static auto q64 = detail::gen_square_table<64>();
  static auto q65 = detail::gen_square_table<65>();
  if (q64[n % 64] == 0)
    return false;
  u32 r = n % 45045;
  if (q63[r % 63] == 0)
    return false;
  if (q65[r % 65] == 0)
    return false;
  if (q11[r % 11] == 0)
    return false;
  u32 sn = sqrt64(n);
  return u64(sn) * sn == n;
}

#endif // ALGO_MATH_SQUARE_TEST
