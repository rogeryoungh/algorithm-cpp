#ifndef ALGO_MATH_POLY_NTT_TWISTED_RADIX_2_BASIC
#define ALGO_MATH_POLY_NTT_TWISTED_RADIX_2_BASIC

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

template <static_modint_concept ModT>
auto &prepare_root_twisted_basic(u32 m) {
  using ValueT = typename ModT::ValueT;
  static constexpr ValueT P = ModT::mod();
  static constexpr ValueT g = 3;
  static constexpr ValueT max_bit = ValueT(1) << std::countr_zero(ModT::mod() - 1);

  static std::vector<ModT> rt{1, 1};
  assert(m <= max_bit);
  while (rt.size() < m) {
    u32 n = rt.size();
    rt.resize(n * 2);
    ModT p = ModT(g).pow((P - 1) / n / 2);
    for (u32 i = n; i < n * 2; i += 2) {
      rt[i] = rt[i / 2], rt[i + 1] = p * rt[i];
    }
  }
  return rt;
}

template <static_modint_concept ModT>
static void ntt_twisted_basic(std::span<ModT> f) { // dif
  i32 n = f.size();
  auto &rt = prepare_root_twisted_basic<ModT>(n);
  for (i32 l = n / 2; l > 0; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = f[i + j + l];
        f[i + j] = x + y;
        f[i + j + l] = ModT::submul(x, y, rt[j + l]);
      }
    }
  }
}

template <static_modint_concept ModT>
static void intt_twisted_basic(std::span<ModT> f) { // dit
  i32 n = f.size();
  auto &rt = prepare_root_twisted_basic<ModT>(n);
  for (i32 l = 1; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = rt[j + l] * f[i + j + l];
        f[i + j] = x + y;
        f[i + j + l] = x - y;
      }
    }
  }
  const ModT ivn = ModT(n).inv();
  for (i32 i = 0; i < n; i++)
    f[i] *= ivn;
  std::reverse(f.begin() + 1, f.end());
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_TWISTED_RADIX_2_BASIC
