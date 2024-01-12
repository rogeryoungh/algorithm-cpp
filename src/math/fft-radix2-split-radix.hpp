#ifndef ALGO_H_MATH_FFT_RADIX2_SPLITRADIX
#define ALGO_H_MATH_FFT_RADIX2_SPLITRADIX

#include "./fft-base.hpp"
#include <vector>
#include <numbers>

ALGO_BEGIN_NAMESPACE

struct FFTRadix2Split : FFTBase {
  std::vector<CP64> rt;
  FFTRadix2Split() : rt(2) {
    rt[0] = rt[1] = CP64{1};
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    for (; n != m; n *= 2) {
      if (n < 32) {
        for (u32 i = n; i != n * 2; ++i) {
          rt[i] = CP64::polar(std::numbers::pi * (i - n) / (n * 2));
        }
      } else {
        CP64 w = CP64::polar(std::numbers::pi / (n * 2));
        for (u32 i = n; i != n * 2; i += 2) {
          rt[i] = rt[i / 2];
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
      CP64 ek = f[0];
      CP64 ok = f[1];
      f[0] = ek + ok;
      f[1] = ek - ok;
    } else {
      for (u32 k = 0; k != n / 4; ++k) {
        CP64 f0 = f[k];
        CP64 f2 = f[k + n / 2];
        CP64 f1 = f[k + n / 4];
        CP64 f3 = f[k + 3 * n / 4];
        CP64 f02s = f0 - f2;
        CP64 f13s = (f1 - f3).mulj();
        f[k] = f0 + f2;
        f[k + n / 2] = CP64::cmul(f02s - f13s, rt[n / 4 + k]);
        f[k + n / 4] = f1 + f3;
        f[k + 3 * n / 4] = (f02s + f13s) * rt[n / 4 + k];
      }
      fft_rec_t<n / 2>(f);
      fft_rec_t<n / 4>(f + n / 2);
      fft_rec_t<n / 4>(f + 3 * n / 4);
    }
  }
  template <u32 n>
  inline void ifft_rec_t(CP64 *f) {
    if constexpr (n <= 1) {
      return;
    } else if constexpr (n == 2) {
      CP64 ek = f[0];
      CP64 ok = f[1];
      f[0] = ek + ok;
      f[1] = ek - ok;
    } else {
      ifft_rec_t<n / 2>(f);
      ifft_rec_t<n / 4>(f + n / 2);
      ifft_rec_t<n / 4>(f + 3 * n / 4);
      for (u32 k = 0; k != n / 4; ++k) {
        CP64 f0 = f[k];
        CP64 f2 = f[k + n / 2] * rt[n / 4 + k];
        CP64 f1 = f[k + n / 4];
        CP64 f3 = CP64::cmul(f[k + 3 * n / 4], rt[n / 4 + k]);
        CP64 f23a = f2 + f3;
        CP64 f23s = (f2 - f3).mulj();
        f[k] = f0 + f23a;
        f[k + n / 2] = f0 - f23a;
        f[k + n / 4] = f1 + f23s;
        f[k + 3 * n / 4] = f1 - f23s;
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
    prepare_root(n / 2);
    fft_rec<1>(f, n);
  }
  void ifft(CP64 *f, u32 n) {
    prepare_root(n / 2);
    ifft_rec<1>(f, n);
  }
};

ALGO_END_NAMESPACE

#endif
