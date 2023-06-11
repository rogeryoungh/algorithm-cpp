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

#include <iostream>

namespace detail {

template <montgomery_modint_concept X8>
struct NttTwistedInfoAvx {
  // using X8 = simd::M32x8<ModT>;
  using ModT = typename X8::ModT;
  using ValueT = typename ModT::ValueT;

  inline static std::vector<X8> rt;
  ValueT P, g, max_bit;

  NttTwistedInfoAvx() {
    P = ModT::mod();
    g = 3;
    max_bit = ValueT(1) << std::countr_zero(ModT::mod() - 1);
    prepare_root<true>(64);
  }

  template <bool init = false>
  void prepare_root(i32 m) {
    if constexpr (init) {
      assert(m <= max_bit);
      rt.resize(m / 8);
      auto rt0 = reinterpret_cast<ModT *>(rt.data());
      rt0[0] = rt0[1] = 1;
      for (i32 n = 2; n < m; n *= 2) {
        ModT p = ModT(g).pow((P - 1) / n / 2);
        for (i32 i = n; i < n * 2; i += 2) {
          rt0[i] = rt0[i / 2], rt0[i + 1] = p * rt0[i];
        }
      }
    } else {
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
  }
  X8 rt_small(i32 L) {
    auto rt0 = reinterpret_cast<ModT *>(rt.data());
    alignas(32) std::array<ModT, 8> r;
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
  using X8 = simd::M32x8<ModT, aligned>;

  // std::cout << "ntt_twisted_avx" << std::endl;

  auto *f = reinterpret_cast<simd::I256 *>(f0.data());
  static auto &info = NttTwistedInfoAvx<X8>::instance();
  static X8 rt2 = info.rt_small(2), rt4 = info.rt_small(4);

  // std::cout << "instance" << std::endl;

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);
  info.prepare_root(n);

  // std::cout << "prepare_root" << std::endl;

  for (i32 l = n / 2; l > 0; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        auto px = &f[i + j], py = &f[i + j + l];
        X8 fx = X8::load(px), fy = X8::load(py);
        X8 rx = fx + fy;
        X8 ry = X8::submul(fx, fy, info.rt[j + l]);
        rx.store(px), ry.store(py);
      }
    }
  }
  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::load(&f[i]);
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi *= rt4;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= rt2;
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi.store(&f[i]);
  }
}

template <class ModT, bool aligned>
static void intt_twisted_avx(std::span<ModT> f0) { // dit
  using X8 = simd::M32x8<ModT, aligned>;

  auto *f = reinterpret_cast<simd::I256 *>(f0.data());
  static auto &info = NttTwistedInfoAvx<X8>::instance();
  static X8 rt2 = info.rt_small(2), rt4 = info.rt_small(4);

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);
  info.prepare_root(n);

  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::load(&f[i]);
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi *= rt2;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= rt4;
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi.store(&f[i]);
  }
  for (i64 l = 1; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        auto px = &f[i + j], py = &f[i + j + l];
        X8 fx = X8::load(px), fy = X8::load(py) * info.rt[j + l];
        X8 rx = fx + fy;
        X8 ry = fx - fy;
        rx.store(px), ry.store(py);
      }
    }
  }
  X8 ivn8 = X8::from(ModT(n8).inv().raw());
  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::load(&f[i]);
    fi *= ivn8;
    fi.store(&f[i]);
  }
  std::reverse(f0.begin() + 1, f0.end());
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_AVX
