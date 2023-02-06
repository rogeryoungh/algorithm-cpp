#ifndef ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_2_BASIC
#define ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_2_BASIC

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

template <static_modint_concept ModT>
struct NttClassicalInfo {
  using ValueT = typename ModT::ValueT;
  static constexpr ValueT P = ModT::mod();
  static constexpr ValueT g = 3;
  static constexpr i32 rank2 = std::countr_zero(P - 1);
  std::array<ModT, rank2 + 1> rt, irt;
  std::array<ModT, std::max<i32>(0, rank2 - 1)> rate2, irate2;

  constexpr NttClassicalInfo() {
    rt[rank2] = ModT(g).pow((P - 1) >> rank2);
    irt[rank2] = rt[rank2].inv();
    for (i32 i = rank2; i >= 1; --i) {
      rt[i - 1] = rt[i] * rt[i];
      irt[i - 1] = irt[i] * irt[i];
    }
    ModT prod = 1, iprod = 1;
    for (i32 i = 0; i < rank2 - 1; ++i) {
      rate2[i] = prod * rt[i + 2];
      irate2[i] = iprod * irt[i + 2];
      prod *= irt[i + 2];
      iprod *= rt[i + 2];
    }
  }
};

template <static_modint_concept ModT>
static void ntt_classical_basic(std::span<ModT> f) { // dif
  static constexpr NttClassicalInfo<ModT> info;
  i32 n = f.size();
  for (i32 l = n / 2; l > 0; l /= 2) {
    ModT r = 1;
    for (i32 i = 0, k = 0; i < n; i += l * 2, ++k) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = f[i + j + l] * r;
        f[i + j] = x + y;
        f[i + j + l] = x - y;
      }
      r *= info.rate2[std::countr_one<u32>(k)];
    }
  }
}

template <static_modint_concept ModT>
static void intt_classical_basic(std::span<ModT> f) { // dit
  static constexpr NttClassicalInfo<ModT> info;
  i32 n = f.size();
  for (i32 l = 1; l < n; l *= 2) {
    ModT r = 1;
    for (i32 i = 0, k = 0; i < n; i += l * 2, ++k) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = f[i + j + l];
        f[i + j] = x + y;
        f[i + j + l] = r * (x - y);
      }
      r *= info.irate2[std::countr_one<u32>(k)];
    }
  }
  const ModT ivn = ModT(n).inv();
  for (i32 i = 0; i < n; i++)
    f[i] *= ivn;
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_2_BASIC
