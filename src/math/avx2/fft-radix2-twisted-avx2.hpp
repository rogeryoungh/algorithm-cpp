#ifndef ALGO_H_MATH_AVX2_FFT_RADIX2_TWISTED
#define ALGO_H_MATH_AVX2_FFT_RADIX2_TWISTED

#include "../../base.hpp"
#include "./complex64x2.hpp"
#include "../../other/align-alloc.hpp"
#include <numbers>

ALGO_BEGIN_NAMESPACE

struct FftR2Avx2T {
  AVec<CP64> rt;
  FftR2Avx2T() : rt(2) {
    rt[0] = rt[1] = CP64{1};
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    for (; n != m; n *= 2) {
      constexpr f64 PI = std::numbers::pi;
      CP64 w = CP64::polar(PI / n);
      for (u32 i = n; i != n * 2; i += 2) {
        rt[i] = rt[i / 2];
        if (i < 32) {
          rt[i + 1] = CP64::polar(PI * (i + 1 - n) / n);
        } else {
          rt[i + 1] = rt[i] * w;
        }
      }
    }
  }
  void fft_small(CP64 *f, u32 n) {
    for (u32 l = n / 2; l != 0; l /= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        for (u32 j = 0; j != l; ++j) {
          CP64 x = f[i + j], y = f[i + j + l];
          f[i + j] = x + y;
          f[i + j + l] = (x - y) * rt[l + j];
        }
      }
    }
  }
  void ifft_small(CP64 *f, u32 n) {
    for (u32 l = 1; l != n; l *= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        for (u32 j = 0; j != l; ++j) {
          CP64 x = f[i + j], y = CP64::cmul(f[i + j + l], rt[l + j]);
          f[i + j] = x + y;
          f[i + j + l] = x - y;
        }
      }
    }
  }
  void fft_butterfly(CP64x2 *f, u32 l, CP64x2 *w) {
    for (u32 j = 0; j != l; ++j) {
      CP64x2 x = f[j], y = f[j + l];
      f[j] = x + y;
      f[j + l] = (x - y) * w[j];
    }
  }
  void fft_layer_last(CP64x2 *f, u32 m) {
    for (u32 j = 0; j != m; ++j) {
      f64x4 abcd = f[j];
      f64x4 baba = _mm256_permute4x64_pd(abcd, 0x44);
      f64x4 dcdc = _mm256_permute4x64_pd(abcd, 0xee);
      f64x4 add = _mm256_add_pd(baba, dcdc);
      f64x4 sub = _mm256_sub_pd(baba, dcdc);
      f[j] = _mm256_blend_pd(add, sub, 0xc);
    }
  }
  void ifft_butterfly(CP64x2 *f, u32 l, CP64x2 *w) {
    for (u32 j = 0; j != l; ++j) {
      CP64x2 x = f[j], y = CP64x2::cmul(f[j + l], w[j]);
      f[j] = x + y;
      f[j + l] = x - y;
    }
  }
  void fft_base(CP64x2 *f, u32 m, CP64x2 *rt2) {
    for (u32 l = m / 2; l != 0; l /= 2) {
      for (u32 i = 0; i != m; i += l * 2) {
        fft_butterfly(f + i, l, rt2 + l);
      }
    }
  }
  void ifft_base(CP64x2 *f, u32 m, CP64x2 *rt2) {
    for (u32 l = 1; l != m; l *= 2) {
      for (u32 i = 0; i != m; i += l * 2) {
        ifft_butterfly(f + i, l, rt2 + l);
      }
    }
  }
  void fft_rec(CP64x2 *f, u32 m, CP64x2 *rt2) {
    constexpr u32 N = 1 << 6;
    if (m <= N) {
      fft_base(f, m, rt2);
    } else {
      u32 l = m / 2;
      fft_butterfly(f, l, rt2 + l);
      fft_rec(f + 0, l, rt2);
      fft_rec(f + l, l, rt2);
    }
  }
  void ifft_rec(CP64x2 *f, u32 n, CP64x2 *rt2) {
    constexpr u32 N = 1 << 6;
    if (n <= N) {
      ifft_base(f, n, rt2);
    } else {
      u32 l = n / 2;
      ifft_rec(f + 0, l, rt2);
      ifft_rec(f + l, l, rt2);
      ifft_butterfly(f, l, rt2 + l);
    }
  }
  void fft(void *p, u32 n) {
    prepare_root(n);
    if (n <= 64)
      return fft_small(reinterpret_cast<CP64 *>(p), n);
    // fft_rec(f, n);
    u32 m = n / 2;
    auto *f = reinterpret_cast<CP64x2 *>(p);
    auto *rt2 = reinterpret_cast<CP64x2 *>(rt.data());
    fft_rec(f, m, rt2);
    fft_layer_last(f, m);
  }
  void ifft(void *p, u32 n) {
    prepare_root(n);
    if (n <= 64)
      return ifft_small(reinterpret_cast<CP64 *>(p), n);
    u32 m = n / 2;
    auto *f = reinterpret_cast<CP64x2 *>(p);
    auto *rt2 = reinterpret_cast<CP64x2 *>(rt.data());
    fft_layer_last(f, m);
    ifft_rec(f, m, rt2);
  }
  void dot(void *p1, const void *p2, u32 n) {
    if (n <= 16) {
      auto *f = reinterpret_cast<CP64 *>(p1);
      auto *g = reinterpret_cast<const CP64 *>(p2);
      for (u32 i = 0; i != n; ++i)
        f[i] *= g[i];
    } else {
      auto *f2 = reinterpret_cast<CP64x2 *>(p1);
      auto *g2 = reinterpret_cast<const CP64x2 *>(p2);
      for (u32 i = 0; i != n / 2; ++i)
        f2[i] *= g2[i];
    }
  }
  void dot2(void *p, u32 n) {
    if (n <= 16) {
      auto *f = reinterpret_cast<CP64 *>(p);
      for (u32 i = 0; i != n; ++i)
        f[i] /= n;
    } else {
      auto *f2 = reinterpret_cast<CP64x2 *>(p);
      for (u32 i = 0; i != n / 2; ++i)
        f2[i] /= n;
    }
  }
};

ALGO_END_NAMESPACE

#endif
