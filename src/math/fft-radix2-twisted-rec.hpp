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
  template <u32 n>
  inline void fft_rec_t(CP64 *f) {
    if constexpr (n <= 1) {
      return;
    } else if constexpr (n == 2) {
      CP64 f0 = f[0];
      CP64 f1 = f[1];
      f[0] = f0 + f1;
      f[1] = f0 - f1;
    } else {
      for (u32 k = 0; k != n / 2; ++k) {
        CP64 f0 = f[k];
        CP64 f1 = f[k + n / 2];
        f[k] = f0 + f1;
        f[k + n / 2] = CP64::cmul(f0 - f1, rt[n / 2 + k]);
      }
      fft_rec_t<n / 2>(f);
      fft_rec_t<n / 2>(f + n / 2);
    }
  }
  template <u32 n>
  inline void ifft_rec_t(CP64 *f) {
    if constexpr (n <= 1) {
      return;
    } else if constexpr (n == 2) {
      CP64 f0 = f[0];
      CP64 f1 = f[1];
      f[0] = f0 + f1;
      f[1] = f0 - f1;
    } else {
      ifft_rec_t<n / 2>(f);
      ifft_rec_t<n / 2>(f + n / 2);
      for (u32 k = 0; k != n / 2; ++k) {
        CP64 f0 = f[k];
        CP64 f1 = f[k + n / 2] * rt[n / 2 + k];
        f[k] = f0 + f1;
        f[k + n / 2] = f0 - f1;
      }
    }
  }
  template <u32 n>
  void fft_rec(CP64 *f, u32 len) {
    if constexpr (n > 0) {
      if (n == len) {
        fft_rec_t<n>(f);
      } else {
        fft_rec<n * 2>(f, len);
      }
    }
  }
  template <u32 n>
  void ifft_rec(CP64 *f, u32 len) {
    if constexpr (n > 0) {
      if (n == len) {
        ifft_rec_t<n>(f);
      } else {
        ifft_rec<n * 2>(f, len);
      }
    }
  }
  void fft(CP64 *f, u32 n) {
    prepare_root(n);
    fft_rec<1>(f, n);
  }
  void ifft(CP64 *f, u32 n) {
    prepare_root(n);
    ifft_rec<1>(f, n);
  }
};

ALGO_END_NAMESPACE

#endif
