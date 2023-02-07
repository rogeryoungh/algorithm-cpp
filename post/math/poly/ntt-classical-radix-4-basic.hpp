#ifndef ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_4_BASIC
#define ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_4_BASIC

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

template <static_modint_concept ModT>
struct NttClassicalInfo4 {
  using ValueT = typename ModT::ValueT;
  static constexpr ValueT P = ModT::mod();
  static constexpr ValueT g = 3;
  static constexpr i32 rank2 = std::countr_zero(P - 1);
  std::array<ModT, rank2 + 1> rt, irt;
  std::array<ModT, std::max<i32>(0, rank2 - 1)> rate2, irate2;
  std::array<ModT, std::max<i32>(0, rank2 - 2)> rate3, irate3;

  constexpr NttClassicalInfo4() {
    rt[rank2] = ModT(g).pow((P - 1) >> rank2);
    irt[rank2] = rt[rank2].inv();
    for (i32 i = rank2; i >= 1; --i) {
      rt[i - 1] = rt[i] * rt[i];
      irt[i - 1] = irt[i] * irt[i];
    }
    ModT prod = 1, iprod = 1;
    for (i32 i = 0; i < rate2.size(); ++i) {
      rate2[i] = prod * rt[i + 2];
      irate2[i] = iprod * irt[i + 2];
      prod *= irt[i + 2];
      iprod *= rt[i + 2];
    }
    prod = 1, iprod = 1;
    for (i32 i = 0; i < rate3.size(); ++i) {
      rate3[i] = prod * rt[i + 3];
      irate3[i] = iprod * irt[i + 3];
      prod *= irt[i + 3];
      iprod *= rt[i + 3];
    }
  }
};

template <static_modint_concept ModT>
static void ntt_classical_basic4(std::span<ModT> f) { // dif
  static constexpr NttClassicalInfo4<ModT> info;
  i32 n = f.size(), l = n / 2, n_4b = std::countr_zero<u32>(n) & 1;
  if (n_4b) {
    for (i32 j = 0; j < l; ++j) {
      ModT x = f[j], y = f[j + l];
      f[j] = x + y;
      f[j + l] = x - y;
    }
    l /= 2;
  }
  for (l /= 2; l >= 1; l /= 4) {
    ModT r = 1, img = info.rt[2];
    for (i32 i = 0, k = 0; i < n; i += l * 4, ++k) {
      ModT r2 = r * r, r3 = r2 * r;
      for (i32 j = 0; j < l; ++j) {
        ModT x0 = f[i + j + 0 * l];
        ModT x1 = f[i + j + 1 * l] * r;
        ModT x2 = f[i + j + 2 * l] * r2;
        ModT x3 = f[i + j + 3 * l] * r3;
        ModT x1x3 = (x1 - x3) * img;
        f[i + j + 0 * l] = x0 + x2 + x1 + x3;
        f[i + j + 1 * l] = x0 + x2 - x1 - x3;
        f[i + j + 2 * l] = x0 - x2 + x1x3;
        f[i + j + 3 * l] = x0 - x2 - x1x3;
      }
      r *= info.rate3[std::countr_one<u32>(k)];
    }
  }
}

template <static_modint_concept ModT>
static void intt_classical_basic4(std::span<ModT> f) { // dit
  static constexpr NttClassicalInfo4<ModT> info;
  i32 n = f.size(), l = 1, n_4b = std::countr_zero<u32>(n) & 1;
  for (; l < (n_4b ? n / 2 : n); l *= 4) {
    ModT r = 1, img = info.irt[2];
    for (i32 i = 0, k = 0; i < n; i += l * 4, ++k) {
      ModT r2 = r * r, r3 = r2 * r;
      for (i32 j = 0; j < l; ++j) {
        ModT x0 = f[i + j + 0 * l];
        ModT x1 = f[i + j + 1 * l];
        ModT x2 = f[i + j + 2 * l];
        ModT x3 = f[i + j + 3 * l];
        ModT x2x3 = (x2 - x3) * img;
        f[i + j + 0 * l] = x0 + x1 + x2 + x3;
        f[i + j + 1 * l] = (x0 - x1 + x2x3) * r;
        f[i + j + 2 * l] = (x0 + x1 - x2 - x3) * r2;
        f[i + j + 3 * l] = (x0 - x1 - x2x3) * r3;
      }
      r *= info.irate3[std::countr_one<u32>(k)];
    }
  }
  if (n_4b) {
    for (i32 j = 0; j < l; ++j) {
      ModT x = f[j], y = f[j + l];
      f[j] = x + y;
      f[j + l] = x - y;
    }
  }
  const ModT ivn = ModT(n).inv();
  for (i32 i = 0; i < n; i++)
    f[i] *= ivn;
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_4_BASIC
