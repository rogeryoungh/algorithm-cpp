#ifndef ALGO_MATH_POLY_NTT_AVX
#define ALGO_MATH_POLY_NTT_AVX

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

#include "../../other/modint/montgomery-x8.hpp"

namespace detail {

template <montgomery_modint_concept ModT>
struct NttTwistedInfoAvx {
  using X8 = simd::M32x8<ModT>;
  using ValueT = typename ModT::ValueT;

  static constexpr ValueT P = ModT::mod();
  static constexpr ValueT g = 3;
  static constexpr ValueT max_bit = ValueT(1) << std::countr_zero(ModT::mod() - 1);

  std::vector<X8> rt;

  NttTwistedInfoAvx() {
    init_rt(64);
  }

  void init_rt(i32 m) {
    assert(m <= max_bit);
    rt.resize(m / 8);
    std::span<ModT> rt0 = as_modt();
    rt0[0] = rt0[1] = 1;
    for (i32 n = 2; n < m; n *= 2) {
      ModT p = ModT(g).pow((P - 1) / n / 2);
      for (i32 i = n; i < n * 2; i += 2) {
        rt0[i] = rt0[i / 2], rt0[i + 1] = p * rt0[i];
      }
    }
  }

  auto as_modt() {
    return std::span{(ModT *)rt.data(), rt.size() * 8};
  }

  void prepare_root(i32 m) {
    assert(m <= max_bit);
    while (rt.size() < m) {
      u32 n = rt.size();
      rt.resize(n * 2);
      ModT p = ModT(g).pow((P - 1) / n / 2 / 8);
      alignas(32) std::array<ModT, 8> arr{};
      for (i32 i = 0; i < 8; ++i)
        arr[i] = i == 0 ? 1 : arr[i - 1] * p;
      X8 pp = arr, p8 = X8::from((arr[7] * p).raw());
      for (i32 i = n; i < n * 2; ++i) {
        rt[i] = pp, pp *= p8;
      }
    }
  }
  template <i32 L>
  X8 rt_small() {
    std::array<ModT, 8> r;
    std::span<ModT> rt0 = {(ModT *)rt.data(), 64};
    std::fill(r.begin(), r.end(), rt0[L + 0]);
    for (i32 i = 0; i < 8; i += L * 2) {
      for (i32 j = 0; j < L; ++j)
        r[i + j + L] = rt0[L + j];
    }
    return r;
  }
  static NttTwistedInfoAvx &instance() {
    static NttTwistedInfoAvx info{};
    return info;
  }
};

template <montgomery_modint_concept ModT, bool aligned>
static void ntt_twisted_avx(std::span<ModT> f0) { // dif
  using X8 = simd::M32x8<ModT>;

  static auto &info = NttTwistedInfoAvx<ModT>::instance();

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);
  std::span<simd::I256> f{(simd::I256 *)f0.data(), u32(n)};
  info.prepare_root(n);

  static X8 rt2 = info.template rt_small<2>();
  static X8 rt4 = info.template rt_small<4>();

  for (i32 l = n / 2; l > 0; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        X8 fx = X8::template load<aligned>(&f[i + j]);
        X8 fy = X8::template load<aligned>(&f[i + j + l]);

        X8 rx = fx + fy;
        X8 ry = info.rt[j + l] * (fx - fy);
        rx.template store<aligned>(&f[i + j]);
        ry.template store<aligned>(&f[i + j + l]);
      }
    }
  }
  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::template load<aligned>(&f[i]);
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi *= rt4;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= rt2;
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi.template store<aligned>(&f[i]);
  }
}

template <montgomery_modint_concept ModT, bool aligned>
static void intt_twisted_avx(std::span<ModT> f0) { // dit
  using X8 = simd::M32x8<ModT>;

  static auto &info = NttTwistedInfoAvx<ModT>::instance();

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);
  std::span<simd::I256> f{(simd::I256 *)f0.data(), u32(n)};
  info.prepare_root(n);

  static X8 rt2 = info.template rt_small<2>();
  static X8 rt4 = info.template rt_small<4>();

  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::template load<aligned>(&f[i]);
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi *= rt2;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= rt4;
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi.template store<aligned>(&f[i]);
  }
  for (i64 l = 1; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        X8 fx = X8::template load<aligned>(&f[i + j]);
        X8 fy = X8::template load<aligned>(&f[i + j + l]) * info.rt[j + l];
        X8 rx = fx + fy;
        X8 ry = fx - fy;
        rx.template store<aligned>(&f[i + j]);
        ry.template store<aligned>(&f[i + j + l]);
      }
    }
  }
  X8 ivn8 = X8::from(ModT(n8).inv().raw());
  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::template load<aligned>(&f[i]);
    fi *= ivn8;
    fi.template store<aligned>(&f[i]);
  }
  std::reverse(f0.begin() + 1, f0.end());
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_AVX
