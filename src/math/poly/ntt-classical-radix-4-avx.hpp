#ifndef ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_4_AVX
#define ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_4_AVX

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <algorithm>
#include <bit>
#include <cassert>
#include <span>
#include <vector>

#include <iostream>

#include "../../other/modint/montgomery-x8.hpp"

namespace detail {

template <montgomery_modint_concept ModT>
struct NttClassicalInfoAvx4 {
  using X8 = simd::M32x8<ModT>;
  using ValueT = typename ModT::ValueT;

  static constexpr ValueT P = ModT::mod();
  static constexpr ValueT g = 3;
  static constexpr i32 rank2 = std::countr_zero(P - 1);
  std::array<ModT, rank2 + 1> rt, irt;
  std::array<ModT, std::max<i32>(0, rank2 - 1)> rate2, irate2;
  std::array<ModT, std::max<i32>(0, rank2 - 2)> rate3, irate3;
  std::array<ModT, std::max<i32>(0, rank2 - 3)> rate4, irate4;
  std::array<X8, std::max<i32>(0, rank2 - 1)> rate2x8, irate2x8;
  std::array<X8, std::max<i32>(0, rank2 - 2)> rate3x8, irate3x8;
  std::array<X8, std::max<i32>(0, rank2 - 3)> rate4ix8, irate4ix8;

  constexpr NttClassicalInfoAvx4() {
    rt[rank2] = ModT(g).pow((P - 1) >> rank2);
    irt[rank2] = rt[rank2].inv();
    for (i32 i = rank2; i >= 1; --i) {
      rt[i - 1] = rt[i] * rt[i];
      irt[i - 1] = irt[i] * irt[i];
    }
    {
      ModT prod = 1, iprod = 1;
      for (i32 i = 0; i < rate2.size(); ++i) {
        rate2[i] = prod * rt[i + 2];
        irate2[i] = iprod * irt[i + 2];
        prod *= irt[i + 2];
        iprod *= rt[i + 2];
        rate2x8[i] = X8::from(rate2[i]);
        irate2x8[i] = X8::from(irate2[i]);
      }
      prod = 1, iprod = 1;
      for (i32 i = 0; i < rate3.size(); ++i) {
        rate3[i] = prod * rt[i + 3];
        irate3[i] = iprod * irt[i + 3];
        prod *= irt[i + 3];
        iprod *= rt[i + 3];
        rate3x8[i] = X8::from(rate3[i]);
        irate3x8[i] = X8::from(irate3[i]);
      }
      prod = 1, iprod = 1;
      for (i32 i = 0; i < rate4.size(); ++i) {
        rate4[i] = prod * rt[i + 4];
        irate4[i] = iprod * irt[i + 4];
        prod *= irt[i + 4];
        iprod *= rt[i + 4];
        std::array<ModT, 8> buf, ibuf;
        for (i32 j = 0; j < 8; ++j) {
          buf[j] = rate4[i].pow(j);
          ibuf[j] = irate4[i].pow(j);
        }
        rate4ix8[i] = buf;
        irate4ix8[i] = ibuf;
      }
    }
  }
  template <i32 L>
  X8 rt_small() {
    std::array<ModT, 8> r;
    std::fill(r.begin(), r.end(), 1);
    if constexpr (L == 2) {
      r[3] = r[7] = rt[2];
    } else if constexpr (L == 4) {
      for (i32 i = 5; i < 8; ++i)
        r[i] = r[i - 1] * rt[3];
    }
    return r;
  }
  template <i32 L>
  X8 irt_small() {
    std::array<ModT, 8> r;
    std::fill(r.begin(), r.end(), 1);
    if constexpr (L == 2) {
      r[3] = r[7] = irt[2];
    } else if constexpr (L == 4) {
      for (i32 i = 5; i < 8; ++i)
        r[i] = r[i - 1] * irt[3];
    }
    return r;
  }
};

template <montgomery_modint_concept ModT, bool aligned>
static void ntt_classical_avx4(std::span<ModT> f0) { // dif
  using X8 = simd::M32x8<ModT, aligned>;

  static NttClassicalInfoAvx4<ModT> info;

  i32 n8 = f0.size(), n = n8 / 8, l = n / 2, n_4b = std::countr_zero<u32>(n) & 1;
  assert(n8 % 16 == 0);
  std::span<simd::I256> f{(simd::I256 *)f0.data(), u32(n)};

  static X8 rt2 = info.template rt_small<2>();
  static X8 rt4 = info.template rt_small<4>();

  if (n_4b) {
    for (i32 j = 0; j < l; ++j) {
      X8 fx = X8::load(&f[j]);
      X8 fy = X8::load(&f[j + l]);
      X8 rx = fx + fy;
      X8 ry = fx - fy;
      rx.store(&f[j]);
      ry.store(&f[j + l]);
    }
    l /= 2;
  }

  for (l /= 2; l >= 1; l /= 4) {
    X8 r = X8::from(ModT(1)), img = X8::from(info.rt[2]);
    for (i32 i = 0, k = 0; i < n; i += l * 4, ++k) {
      X8 r2 = r * r, r3 = r2 * r;
      for (i32 j = 0; j < l; ++j) {
        X8 x0 = X8::load(&f[i + j + 0 * l]);
        X8 x1 = X8::load(&f[i + j + 1 * l]) * r;
        X8 x2 = X8::load(&f[i + j + 2 * l]) * r2;
        X8 x3 = X8::load(&f[i + j + 3 * l]) * r3;
        X8 x1x3 = X8::submul(x1, x3, img);
        X8 x02 = x0 + x2, x0_2 = x0 - x2;
        X8 x13 = x1 + x3;
        X8 y0 = x02 + x13;
        X8 y1 = x02 - x13;
        X8 y2 = x0_2 + x1x3;
        X8 y3 = x0_2 - x1x3;
        y0.store(&f[i + j + 0 * l]);
        y1.store(&f[i + j + 1 * l]);
        y2.store(&f[i + j + 2 * l]);
        y3.store(&f[i + j + 3 * l]);
      }
      r *= info.rate3x8[std::countr_one<u32>(k)];
    }
  }
  X8 rti = X8::from(ModT(1));
  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::load(&f[i]);
    fi *= rti;
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi *= rt4;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= rt2;
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi.store(&f[i]);
    rti *= info.rate4ix8[std::countr_one<u32>(i)];
  }
}

template <montgomery_modint_concept ModT, bool aligned>
static void intt_classical_avx4(std::span<ModT> f0) { // dit
  using X8 = simd::M32x8<ModT, aligned>;

  static NttClassicalInfoAvx4<ModT> info;

  i32 n8 = f0.size(), n = n8 / 8, l = 1, n_4b = std::countr_zero<u32>(n) & 1;
  assert(n8 % 16 == 0);
  std::span<simd::I256> f{(simd::I256 *)f0.data(), u32(n)};

  static X8 rt2 = info.template irt_small<2>();
  static X8 rt4 = info.template irt_small<4>();

  X8 rti = X8::from(ModT(1));
  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::load(&f[i]);
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi *= rt2;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= rt4;
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi *= rti;
    fi.store(&f[i]);
    rti *= info.irate4ix8[std::countr_one<u32>(i)];
  }
  for (; l < (n_4b ? n / 2 : n); l *= 4) {
    X8 r = X8::from(ModT(1)), img = X8::from(info.irt[2]);
    for (i32 i = 0, k = 0; i < n; i += l * 4, ++k) {
      X8 r2 = r * r, r3 = r2 * r;
      for (i32 j = 0; j < l; ++j) {
        X8 x0 = X8::load(&f[i + j + 0 * l]);
        X8 x1 = X8::load(&f[i + j + 1 * l]);
        X8 x2 = X8::load(&f[i + j + 2 * l]);
        X8 x3 = X8::load(&f[i + j + 3 * l]);
        X8 x2x3 = X8::submul(x2, x3, img);
        X8 x01 = x0 + x1, x0_1 = x0 - x1;
        X8 x23 = x2 + x3;
        X8 y0 = x01 + x23;
        X8 y1 = X8::addmul(x0_1, x2x3, r);
        X8 y2 = X8::submul(x01, x23, r2);
        X8 y3 = X8::submul(x0_1, x2x3, r3);
        y0.store(&f[i + j + 0 * l]);
        y1.store(&f[i + j + 1 * l]);
        y2.store(&f[i + j + 2 * l]);
        y3.store(&f[i + j + 3 * l]);
      }
      r *= info.irate3x8[std::countr_one<u32>(k)];
    }
  }
  if (n_4b) {
    for (i32 j = 0; j < l; ++j) {
      X8 fx = X8::load(&f[j]);
      X8 fy = X8::load(&f[j + l]);
      X8 rx = fx + fy;
      X8 ry = fx - fy;
      rx.store(&f[j]);
      ry.store(&f[j + l]);
    }
  }
  X8 ivn8 = X8::from(ModT(n8).inv());
  for (i32 i = 0; i < n; ++i) {
    X8 fi = X8::load(&f[i]);
    fi *= ivn8;
    fi.store(&f[i]);
  }
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_4_AVX
