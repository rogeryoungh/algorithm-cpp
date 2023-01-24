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
  constexpr value P = ModT::get_mod();
  constexpr value g = 3;
  constexpr u64 max_bit = (P - 1) & -(P - 1);

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
};

template <static_modint_concept ModT>
constexpr inline ModT gen_ivn(u32 n) {
  enum : typename ModT::value_type {
    P = ModT::get_mod(),
  };
  return P - (P - 1) / n;
}

template <static_modint_concept ModT>
static void ntt(std::span<ModT> f) { // dif
  u32 n = f.size();
  assert(std::has_single_bit(n));
  auto &rt = prepare_root<ModT>(n);
  detail::ntt_size += n;
  for (u32 l = n / 2; l > 0; l /= 2) {
    for (u32 i = 0; i < n; i += l * 2) {
      for (u32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = f[i + j + l];
        f[i + j + l] = (x - y) * rt[j + l];
        f[i + j] = x + y;
      }
    }
  }
}

template <static_modint_concept ModT>
static void intt(std::span<ModT> f) { // dit
  u32 n = f.size();
  assert(std::has_single_bit(n));
  auto &rt = prepare_root<ModT>(n);
  detail::ntt_size += n;
  for (u32 l = 1; l < n; l *= 2) {
    for (u32 i = 0; i < n; i += l * 2) {
      for (u32 j = 0; j < l; ++j) {
        ModT x = f[i + j], y = rt[j + l] * f[i + j + l];
        f[i + j] = x + y, f[i + j + l] = x - y;
      }
    }
  }
  const ModT ivn = gen_ivn<ModT>(n);
  for (u32 i = 0; i < n; i++)
    f[i] *= ivn;
  std::reverse(f.begin() + 1, f.end());
}

// 这个性能差距太神秘了。。本地测出来能快 10% 左右，OJ 上无区别
template <static_raw32_modint_concept ModT>
static void ntt(std::span<ModT> f_) { // dif
  auto f = std::bit_cast<std::span<u32>>(f_);
  constexpr u32 P = ModT::get_mod();

  u32 n = f.size();
  assert(std::has_single_bit(n));
  auto &rt = prepare_root<ModT>(n);
  detail::ntt_size += n;
  for (u32 l = n / 2; l > 0; l /= 2) {
    for (u32 i = 0; i < n; i += l * 2) {
      for (u32 j = 0; j < l; ++j) {
        u32 x = f[i + j], y = f[i + j + l];
        f[i + j + l] = u64(x - y + P) * rt[j + l].val() % P;
        f[i + j] = detail::mo<P>(x + y);
      }
    }
  }
}

template <static_raw32_modint_concept ModT>
static void intt(std::span<ModT> f_) { // dit
  auto f = std::bit_cast<std::span<u32>>(f_);
  constexpr u32 P = ModT::get_mod();
  constexpr u32 P2 = P * 2;

  u32 n = f.size();
  assert(std::has_single_bit(n));
  auto &rt = prepare_root<ModT>(n);
  detail::ntt_size += n;
  for (u32 l = 1; l < n; l *= 2) {
    for (u32 i = 0; i < n; i += l * 2) {
      for (u32 j = 0; j < l; ++j) {
        u32 x = detail::mo<P2>(f[i + j]), y = u64(rt[j + l].val()) * f[i + j + l] % P;
        f[i + j] = x + y, f[i + j + l] = x - y + P;
      }
    }
  }
  const u32 ivn = P - (P - 1) / n;
  for (u32 i = 0; i < n; i++)
    f[i] = u64(ivn) * f[i] % P;
  std::reverse(f.begin() + 1, f.end());
}

#endif // ALGO_MATH_POLY_NTT
