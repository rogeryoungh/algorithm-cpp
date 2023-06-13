#ifndef ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_2_AVX
#define ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_2_AVX

#include "../../base.hpp"
#include "../../other/modint/montgomery-x8.hpp"
#include "vec-dots.hpp"

#include <bit>
#include <cassert>
#include <span>
#include <vector>

namespace detail {

template <class X8>
struct NttClassicalInfoAvx {
  using ModT = typename X8::ModT;
  using ValueT = typename ModT::ValueT;
  std::array<ModT, 64> rt, irt, rate2, irate2, rate4, irate4;
  std::array<X8, 64> rate4ix8, irate4ix8;
  X8 rt2, irt2, rt4, irt4;

  NttClassicalInfoAvx() {
    const ValueT P = ModT::mod();
    const ValueT g = 3;
    const i32 rank2 = std::countr_zero(P - 1);
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
    prod = 1, iprod = 1;
    alignas(32) std::array<ModT, 8> r, ir;
    for (i32 i = 0; i < rank2 - 3; ++i) {
      rate4[i] = prod * rt[i + 4];
      irate4[i] = iprod * irt[i + 4];
      prod *= irt[i + 4];
      iprod *= rt[i + 4];
      for (i32 j = 0; j < 8; ++j) {
        r[j] = rate4[i].pow(j);
        ir[j] = irate4[i].pow(j);
      }
      rate4ix8[i] = r;
      irate4ix8[i] = ir;
    }
    std::fill(r.begin(), r.end(), 1);
    std::fill(ir.begin(), ir.end(), 1);
    r[3] = r[7] = rt[2];
    ir[3] = ir[7] = irt[2];
    rt2 = r, irt2 = ir;
    std::fill(r.begin(), r.end(), 1);
    std::fill(ir.begin(), ir.end(), 1);
    for (i32 i = 5; i < 8; ++i) {
      r[i] = r[i - 1] * rt[3];
      ir[i] = ir[i - 1] * irt[3];
    }
    rt4 = r, irt4 = ir;
  }
};

template <montgomery_modint_concept ModT, bool aligned>
static void ntt_classical_avx(std::span<ModT> f0) { // dif
  using X8 = simd::M32x8<ModT>;
  auto *f = std::bit_cast<X8 *>(f0.data());
  static const NttClassicalInfoAvx<X8> info;

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);

  for (i32 l = n / 2; l >= 1 * 1; l /= 2) {
    ModT r = 1;
    for (i32 i = 0, k = 0; i < n; i += l * 2, ++k) {
      X8 r8 = X8::from(r);
      for (i32 j = 0; j < l; ++j) {
        X8 fx = f[i + j];
        X8 fy = f[i + j + l] * r8;
        f[i + j] = fx + fy;
        f[i + j + l] = fx - fy;
      }
      r *= info.rate2[std::countr_one<u32>(k)];
    }
  }
  X8 rti = X8::from(ModT(1));
  for (i32 i = 0; i < n; ++i) {
    X8 fi = f[i];
    fi *= rti;
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi *= info.rt4;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= info.rt2;
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    f[i] = fi;
    rti *= info.rate4ix8[std::countr_one<u32>(i)];
  }
}

template <montgomery_modint_concept ModT, bool aligned>
static void intt_classical_avx(std::span<ModT> f0) { // dit
  using X8 = simd::M32x8<ModT>;
  auto *f = std::bit_cast<X8 *>(f0.data());
  static const NttClassicalInfoAvx<X8> info;

  i32 n8 = f0.size(), n = n8 / 8;
  assert(n8 % 16 == 0);

  X8 rti = X8::from(ModT(1));
  for (i32 i = 0; i < n; ++i) {
    X8 fi = f[i];
    fi = fi.template neg<0b10101010>() + fi.template shuffle<0b10110001>();
    fi *= info.irt2;
    fi = fi.template neg<0b11001100>() + fi.template shuffle<0b01001110>();
    fi *= info.irt4;
    fi = fi.template neg<0b11110000>() + fi.template shufflex4<0b01>();
    fi *= rti;
    f[i] = fi;
    rti *= info.irate4ix8[std::countr_one<u32>(i)];
  }
  for (i64 l = 1; l < n; l *= 2) {
    ModT r = 1;
    for (i32 i = 0, k = 0; i < n; i += l * 2, ++k) {
      X8 r8 = X8::from(r);
      for (i32 j = 0; j < l; ++j) {
        X8 fx = f[i + j];
        X8 fy = f[i + j + l];
        f[i + j] = fx + fy;
        f[i + j + l] = X8::submul(fx, fy, r8);
      }
      r *= info.irate2[std::countr_one<u32>(k)];
    }
  }
  dot_v<ModT, aligned, false>(f0, ModT(n8).inv());
}

} // namespace detail

#endif // ALGO_MATH_POLY_NTT_CLASSICAL_RADIX_2_AVX
