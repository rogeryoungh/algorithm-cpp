#ifndef ALGO_MATH_POLY_NTT
#define ALGO_MATH_POLY_NTT

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

u32 ntt_size = 0;

template <u32 M>
inline u32 mo(u32 x) {
  return x >= M ? x - M : x;
}

} // namespace detail

template <static_modint_concept ModT>
auto &prepare_root(u32 m) {
  using value = typename ModT::value_type;
  static constexpr value P = ModT::get_mod();
  static constexpr value g = 3;
  static constexpr u64 max_bit = (P - 1) & -(P - 1);

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
static void ntt(std::span<ModT> f) { // dif
  i32 n = f.size();
  assert(std::has_single_bit<u32>(n));
  auto &rt = prepare_root<ModT>(n);
  detail::ntt_size += n;
  for (i32 l = n / 2; l > 0; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = f[i + j + l];
        f[i + j + l] = rt[j + l] * (x - y);
        f[i + j] = x + y;
      }
    }
  }
}

template <static_modint_concept ModT>
static void intt(std::span<ModT> f) { // dit
  i32 n = f.size();
  assert(std::has_single_bit<u32>(n));
  auto &rt = prepare_root<ModT>(n);
  detail::ntt_size += n;
  for (i32 l = 1; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = rt[j + l] * f[i + j + l];
        f[i + j] = x + y, f[i + j + l] = x - y;
      }
    }
  }
  const ModT ivn = ModT::get_mod() - (ModT::get_mod() - 1) / n;
  for (i32 i = 0; i < n; i++)
    f[i] *= ivn;
  std::reverse(f.begin() + 1, f.end());
}

template <u32 MOD>
struct BarrettMul {
  u32 b = 0, c = 0;
  BarrettMul() = default;
  inline BarrettMul(u32 x) : b(x), c((u64(x) << 32) / MOD) {}
  inline u32 operator*(u32 a) { // <= 2M !!!
    return u64(a) * b - (u64(a) * c >> 32) * MOD;
  }
};

template <static_modint_concept ModT>
auto &prepare_root_barrett(u32 m) {
  using value = typename ModT::value_type;
  static constexpr value P = ModT::get_mod();
  static constexpr value g = 3;
  static constexpr u64 max_bit = (P - 1) & -(P - 1);

  static std::vector<BarrettMul<P>> rt{{1}, {1}};
  assert(m <= max_bit);
  while (rt.size() < m) {
    u32 n = rt.size();
    rt.resize(n * 2);
    u32 p = ModT(g).pow((P - 1) / n / 2).val();
    for (u32 i = n; i < n * 2; i += 2) {
      rt[i] = rt[i / 2], rt[i + 1] = u64(p) * rt[i].b % P;
    }
  }
  return rt;
}

template <static_raw32_modint_concept ModT>
static void ntt(std::span<ModT> f_) { // dif
  i32 n = f_.size();
  assert(std::has_single_bit<u32>(n));
  auto &rt = prepare_root_barrett<ModT>(n);
  detail::ntt_size += n;

  constexpr u32 P = ModT::get_mod();
  constexpr u32 P2 = P * 2;

  auto f = std::bit_cast<std::span<u32>>(f_);
  for (i32 l = n / 2; l > 0; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        u32 x = f[i + j], y = f[i + j + l];
        f[i + j + l] = rt[j + l] * (x - y + P2);
        f[i + j] = detail::mo<P2>(x + y);
      }
    }
  }
  for (i32 i = 0; i < n; i++)
    f[i] = detail::mo<P>(detail::mo<P2>(f[i]));
}

template <static_raw32_modint_concept ModT>
static void intt(std::span<ModT> f_) { // dit
  i32 n = f_.size();
  assert(std::has_single_bit<u32>(n));
  auto &rt = prepare_root_barrett<ModT>(n);
  detail::ntt_size += n;

  constexpr u32 P = ModT::get_mod();
  constexpr u32 P2 = P * 2;

  auto f = std::bit_cast<std::span<u32>>(f_);
  for (i32 l = 1; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        u32 x = detail::mo<P2>(f[i + j]);
        u32 y = rt[j + l] * f[i + j + l];
        f[i + j] = x + y, f[i + j + l] = x - y + P2;
      }
    }
  }
  const u32 ivn = P - (P - 1) / n;
  for (i32 i = 0; i < n; i++)
    f[i] = f[i] * u64(ivn) % P;
  std::reverse(f.begin() + 1, f.end());
}

#endif // ALGO_MATH_POLY_NTT
