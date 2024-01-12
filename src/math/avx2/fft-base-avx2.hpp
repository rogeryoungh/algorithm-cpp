#ifndef ALGO_H_MATH_AVX2_FFT_BASE
#define ALGO_H_MATH_AVX2_FFT_BASE

#include "../../number/avx2/complex64x2.hpp"

ALGO_BEGIN_NAMESPACE

struct FFT64BaseAVX2 {
  void dot(CP64 *f, const CP64 *g, u32 n) {
    if (n <= 1) {
      for (u32 i = 0; i != n; ++i)
        f[i] *= g[i];
    } else {
      auto *fx = reinterpret_cast<CP64x2 *>(f);
      auto *gx = reinterpret_cast<const CP64x2 *>(g);
      for (u32 i = 0; i != n / 2; ++i)
        fx[i] *= gx[i];
    }
  }
  void rescale(CP64 *f, u32 n) {
    f64 ivn = f64(1) / n;
    if (n <= 1) {
      for (u32 i = 0; i != n; ++i)
        f[i] *= ivn;
    } else {
      auto *fx = reinterpret_cast<CP64x2 *>(f);
      for (u32 i = 0; i != n / 2; ++i)
        fx[i] *= ivn;
    }
  }
};

ALGO_END_NAMESPACE

#endif
