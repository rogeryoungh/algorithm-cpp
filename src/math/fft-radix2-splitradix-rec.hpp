#ifndef ALGO_H_MATH_FFT_RADIX2_TWISTED
#define ALGO_H_MATH_FFT_RADIX2_TWISTED

#include "../base.hpp"
#include "./complex64.hpp"
#include <vector>
#include <numbers>

ALGO_BEGIN_NAMESPACE

struct FFTRadix2Split {
  std::vector<CP64> rt1, rt2;
  FFTRadix2Split() : rt1(2), rt2(2) {
    rt2[0] = rt2[1] = rt1[0] = rt1[1] = CP64{1};
  }
  void prepare_root(u32 m) {
    u32 n = rt1.size();
    if (n >= m)
      return;
    rt1.resize(m);
    rt2.resize(m);
    for (; n != m; n *= 2) {
      if (n < 32) {
        for (u32 i = n; i != n * 2; ++i) {
          rt1[i] = CP64::polar(std::numbers::pi * (i - n) / (n * 2));
          rt2[i] = CP64::polar(std::numbers::pi * (i - n) * 3 / (n * 2));
        }
      } else {
        CP64 w = CP64::polar(std::numbers::pi / (n * 2));
        CP64 w3 = CP64::polar(std::numbers::pi * 3 / (n * 2));
        for (u32 i = n; i != n * 2; i += 2) {
          rt1[i] = rt1[i / 2];
          rt1[i + 1] = rt1[i] * w;
          rt2[i] = rt2[i / 2];
          rt2[i + 1] = rt2[i] * w3;
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
        f[k + n / 2] = CP64::cmul(f02s - f13s, rt1[n / 4 + k]);
        f[k + n / 4] = f1 + f3;
        f[k + 3 * n / 4] = CP64::cmul(f02s + f13s, rt2[n / 4 + k]);
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
        CP64 f2 = f[k + n / 2] * rt1[n / 4 + k];
        CP64 f1 = f[k + n / 4];
        CP64 f3 = f[k + 3 * n / 4] * rt2[n / 4 + k];
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
  void fft_rec(CP64 *f, int len) {
    if constexpr (n <= 1) {
      return;
    } else if (n == len) {
      return fft_rec_t<n>(f);
    } else {
      return fft_rec<n / 2>(f, len);
    }
  }
  template <u32 n>
  void ifft_rec(CP64 *f, int len) {
    if constexpr (n <= 1) {
      return;
    } else if (n == len) {
      return ifft_rec_t<n>(f);
    } else {
      return ifft_rec<n / 2>(f, len);
    }
  }
  void fft(CP64 *f, u32 n) {
    prepare_root(n);
    fft_rec<(1u << 25)>(f, n);
  }
  void ifft(CP64 *f, u32 n) {
    prepare_root(n);
    ifft_rec<(1u << 25)>(f, n);
  }
  void dot(CP64 *f, const CP64 *g, u32 n) {
    for (u32 i = 0; i != n; ++i)
      f[i] *= g[i];
  }
  void div2n(CP64 *f, u32 n) {
    f64 ivn = f64(1) / n;
    for (u32 i = 0; i != n; ++i)
      f[i] *= ivn;
  }
};

ALGO_END_NAMESPACE

#endif
