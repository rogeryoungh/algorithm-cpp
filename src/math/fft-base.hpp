#ifndef ALGO_H_MATH_FFT_BASE
#define ALGO_H_MATH_FFT_BASE

#include "../base.hpp"
#include "../number/complex64.hpp"

ALGO_BEGIN_NAMESPACE

struct FFTBase {
  void dot(CP64 *f, const CP64 *g, u32 n) {
    for (u32 i = 0; i != n; ++i)
      f[i] *= g[i];
  }
  void rescale(CP64 *f, u32 n) {
    f64 ivn = f64(1) / n;
    for (u32 i = 0; i != n; ++i)
      f[i] *= ivn;
  }
};

ALGO_END_NAMESPACE

#endif
