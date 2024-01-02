#ifndef ALGO_H_MATH_FFT_RADIX2
#define ALGO_H_MATH_FFT_RADIX2

#include "../base.hpp"
#include "./complex64.hpp"
#include <array>
#include <numbers>

ALGO_BEGIN_NAMESPACE

struct FftR2 {
  std::array<CP64, 64> rate2;
  FftR2() {
    constexpr u32 rank2 = 22;
    std::array<CP64, rank2 + 1> rt;
    for (u32 i = 0; i != rank2 + 1; ++i) {
      rt[i] = i == 0 ? CP64{1} : CP64::polar(2 * std::numbers::pi / (u64(1) << i));
    }
    CP64 prod = CP64{1};
    for (u32 i = 0; i != rank2 - 1; ++i) {
      rate2[i] = prod * rt[i + 2];
      prod = CP64::cmul(prod, rt[i + 2]);
    }
  }
  void fft(CP64 *f, u32 n) const {
    for (u32 l = n / 2; l != 0; l /= 2) {
      CP64 r = CP64{1};
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          CP64 x = f[i + j], y = f[i + j + l] * r;
          f[i + j] = x + y;
          f[i + j + l] = x - y;
        }
        r *= rate2[std::countr_one(k)];
      }
    }
  }
  void ifft(CP64 *f, u32 n) const {
    for (u32 l = 1; l != n; l *= 2) {
      CP64 r = CP64{1};
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          CP64 x = f[i + j], y = f[i + j + l];
          f[i + j] = x + y;
          f[i + j + l] = CP64::cmul(x - y, r);
        }
        r *= rate2[std::countr_one(k)];
      }
    }
  }
  void dot(CP64 *f, const CP64 *g, u32 n) {
    for (u32 i = 0; i != n; ++i)
      f[i] *= g[i];
  }
  void dot2(CP64 *f, u32 n) {
    for (u32 i = 0; i != n; ++i)
      f[i] /= n;
  }
};

ALGO_END_NAMESPACE

#endif
