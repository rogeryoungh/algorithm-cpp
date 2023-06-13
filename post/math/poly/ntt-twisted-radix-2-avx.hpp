#ifndef ALGO_MATH_POLY_NTT_AVX
#define ALGO_MATH_POLY_NTT_AVX

#include "../../base.hpp"
#include "../../other/modint/montgomery-x8.hpp"
#include "vec-dots.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

template <class X8>
struct NttTwistedInfoAvx {
  using ModT = typename X8::ModT;
  using ValueT = typename ModT::ValueT;

  std::vector<X8> rt;
  ValueT P, g, max_bit;
  X8 rt2, irt2, rt4, irt4;

  NttTwistedInfoAvx() {
    P = ModT::mod();
    g = 3;
    max_bit = ValueT(1) << std::countr_zero(ModT::mod() - 1);
    prepare_root(16);

    auto rt0 = reinterpret_cast<ModT *>(rt.data());

    alignas(32) std::array<ModT, 8> r;

    std::fill(r.begin(), r.end(), rt0[2]);
    r[3] = r[7] = rt0[3];
    rt2 = r;
    std::fill(r.begin(), r.end(), rt0[4]);
    r[5] = rt0[5], r[6] = rt0[6], r[7] = rt0[7];
    rt4 = r;
  }

  void prepare_root(u32 m) {
    m = std::bit_ceil(m);
    assert(m <= max_bit);

    u32 n = rt.size() * 8;
    if (n < m) {
      rt.resize(m / 8);
      auto rt0 = reinterpret_cast<ModT *>(rt.data());
      if (n == 0) {
        n = 2, rt0[0] = rt0[1] = 1;
      }
      for (; n < m; n *= 2) {
        ModT p = ModT(g).pow((P - 1) / n / 2);
        for (u32 i = n; i < n * 2; i += 2) {
          rt0[i] = rt0[i / 2], rt0[i + 1] = p * rt0[i];
        }
      }
    }
  }
  static NttTwistedInfoAvx &instance() {
    static NttTwistedInfoAvx info{};
    return info;
  }
};

template <montgomery_modint_concept ModT, bool aligned>
static void ntt_twisted_avx(std::span<ModT> f0) { // dif
  using X8 = simd::M32x8<ModT>;
  auto *f = std::bit_cast<X8 *>(f0.data());
  static auto &info = NttTwistedInfoAvx<X8>::instance();

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);
  info.prepare_root(n8);

  for (i32 l = n / 2; l > 0; l /= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        X8 fx = f[i + j], fy = f[i + j + l];
        f[i + j] = fx + fy;
        f[i + j + l] = (fx - fy) * info.rt[j + l];
      }
    }
  }
  for (i32 i = 0; i < n; ++i) {
    X8 fi = f[i];
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi *= info.rt4;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= info.rt2;
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    f[i] = fi;
  }
}

template <class ModT, bool aligned>
static void intt_twisted_avx(std::span<ModT> f0) { // dit
  using X8 = simd::M32x8<ModT>;
  auto *f = std::bit_cast<X8 *>(f0.data());
  static auto &info = NttTwistedInfoAvx<simd::M32x8<ModT>>::instance();

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);
  info.prepare_root(n8);

  for (i32 i = 0; i < n; ++i) {
    X8 fi = f[i];
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi *= info.rt2;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= info.rt4;
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    f[i] = fi;
  }
  for (i64 l = 1; l < n; l *= 2) {
    for (i32 i = 0; i < n; i += l * 2) {
      for (i32 j = 0; j < l; ++j) {
        X8 fx = f[i + j], fy = f[i + j + l] * info.rt[j + l];
        f[i + j] = fx + fy;
        f[i + j + l] = fx - fy;
      }
    }
  }
  dot_v<ModT, aligned, false>(f0, ModT(n8).inv());
  std::reverse(f0.begin() + 1, f0.end());
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_AVX
