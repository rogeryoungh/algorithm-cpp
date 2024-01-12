#ifndef ALGO_H_MATH_FFT_RADIX2
#define ALGO_H_MATH_FFT_RADIX2

#include "./fft-base.hpp"
#include <numbers>
#include <vector>

ALGO_BEGIN_NAMESPACE

struct FFTRadix2 : FFTBase {
  std::vector<CP64> rt;
  FFTRadix2() {
    rt = {CP64{1}};
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    for (; n != m; n *= 2) {
      CP64 w = CP64::polar(std::numbers::pi / n / 2);
      for (u32 i = n; i != n * 2; ++i) {
        rt[i] = rt[i - n] * w;
      }
    }
  }
  void fft(CP64 *f, u32 n) {
    prepare_root(n);
    for (u32 l = n / 2; l != 0; l /= 2) {
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          CP64 x = f[i + j], y = f[i + j + l] * rt[k];
          f[i + j] = x + y;
          f[i + j + l] = x - y;
        }
      }
    }
  }
  void ifft(CP64 *f, u32 n) {
    prepare_root(n);
    for (u32 l = 1; l != n; l *= 2) {
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          CP64 x = f[i + j], y = f[i + j + l];
          f[i + j] = x + y;
          f[i + j + l] = CP64::cmul(x - y, rt[k]);
        }
      }
    }
  }
};

ALGO_END_NAMESPACE

#endif
