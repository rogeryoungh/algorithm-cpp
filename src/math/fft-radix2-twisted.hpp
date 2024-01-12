#ifndef ALGO_H_MATH_FFT_RADIX2_TWISTED
#define ALGO_H_MATH_FFT_RADIX2_TWISTED

#include "./fft-base.hpp"
#include <vector>
#include <numbers>

ALGO_BEGIN_NAMESPACE

struct FFTRadix2Twisted : FFTBase {
  std::vector<CP64> rt;
  FFTRadix2Twisted() : rt(2) {
    rt[0] = rt[1] = CP64{1};
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    for (; n != m; n *= 2) {
      CP64 w = CP64::polar(std::numbers::pi / n);
      for (u32 i = n; i != n * 2; i += 2) {
        rt[i] = rt[i / 2];
        if (i < (1 << 4)) {
          rt[i + 1] = CP64::polar(std::numbers::pi * (i + 1 - n) / n);
        } else {
          rt[i + 1] = rt[i] * w;
        }
      }
    }
  }
  void fft_butterfly(CP64 *f, u32 l, CP64 *w) {
    for (u32 j = 0; j != l; ++j) {
      CP64 x = f[j], y = f[j + l];
      f[j] = x + y;
      f[j + l] = CP64::cmul(x - y, w[j]);
    }
  }
  void ifft_butterfly(CP64 *f, u32 l, CP64 *w) {
    for (u32 j = 0; j != l; ++j) {
      CP64 x = f[j], y = f[j + l] * w[j];
      f[j] = x + y;
      f[j + l] = x - y;
    }
  }
  void fft_rec(CP64 *f, u32 n) {
    constexpr u32 N = 1 << 10;
    if (n <= N) {
      for (u32 l = n / 2; l != 0; l /= 2) {
        for (u32 i = 0; i != n; i += l * 2) {
          fft_butterfly(f + i, l, rt.data() + l);
        }
      }
    } else {
      u32 l = n / 2;
      fft_butterfly(f, l, rt.data() + l);
      fft_rec(f + 0, l);
      fft_rec(f + l, l);
    }
  }
  void ifft_rec(CP64 *f, u32 n) {
    constexpr u32 N = 1 << 10;
    if (n <= N) {
      for (u32 l = 1; l != n; l *= 2) {
        for (u32 i = 0; i != n; i += l * 2) {
          ifft_butterfly(f + i, l, rt.data() + l);
        }
      }
    } else {
      u32 l = n / 2;
      ifft_rec(f + 0, l);
      ifft_rec(f + l, l);
      ifft_butterfly(f, l, rt.data() + l);
    }
  }
  void fft(CP64 *f, u32 n) {
    prepare_root(n);
    fft_rec(f, n);
  }
  void ifft(CP64 *f, u32 n) {
    prepare_root(n);
    ifft_rec(f, n);
  }
};

ALGO_END_NAMESPACE

#endif
