#ifndef ALGO_H_MATH_AVX2_FFT_RADIX2_TWISTED
#define ALGO_H_MATH_AVX2_FFT_RADIX2_TWISTED

#include "../../base.hpp"
#include "./complex64x2.hpp"
#include "../../other/align-alloc.hpp"
#include <numbers>

ALGO_BEGIN_NAMESPACE

struct FFT64Radix2TwistedAVX2 {
  AVec<CP64> rt;
  FFT64Radix2TwistedAVX2() : rt(2) {
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
  void fft_rec(CP64x2 *f, u32 m, CP64x2 *rt2) {
    constexpr u32 N = 1 << 6;
    if (m <= N) {
      for (u32 l = m / 2; l != 0; l /= 2) {
        for (u32 i = 0; i != m; i += l * 2) {
          fft_butterfly(f + i, l, rt2 + l);
        }
      }
    } else {
      u32 l = m / 2;
      fft_butterfly(f, l, rt2 + l);
      fft_rec(f + 0, l, rt2);
      fft_rec(f + l, l, rt2);
    }
  }
  void ifft_rec(CP64x2 *f, u32 m, CP64x2 *rt2) {
    constexpr u32 N = 1 << 6;
    if (m <= N) {
      for (u32 l = 1; l != m; l *= 2) {
        for (u32 i = 0; i != m; i += l * 2) {
          ifft_butterfly(f + i, l, rt2 + l);
        }
      }
    } else {
      u32 l = m / 2;
      ifft_rec(f + 0, l, rt2);
      ifft_rec(f + l, l, rt2);
      ifft_butterfly(f, l, rt2 + l);
    }
  }
  void fft(CP64 *f, u32 n) {
    prepare_root(n);
    if (n <= 8)
      return fft_small(f, n);
    u32 m = n / 2;
    auto *fx = reinterpret_cast<CP64x2 *>(f);
    auto *rtx = reinterpret_cast<CP64x2 *>(rt.data());
    fft_rec(fx, m, rtx);
    fft_layer_last(fx, m);
  }
  void ifft(CP64 *f, u32 n) {
    prepare_root(n);
    if (n <= 8)
      return ifft_small(f, n);
    u32 m = n / 2;
    auto *fx = reinterpret_cast<CP64x2 *>(f);
    auto *rtx = reinterpret_cast<CP64x2 *>(rt.data());
    fft_layer_last(fx, m);
    ifft_rec(fx, m, rtx);
  }
  void dot(CP64 *f, const CP64 *g, u32 n) {
    if (n <= 16) {
      for (u32 i = 0; i != n; ++i)
        f[i] *= g[i];
    } else {
      auto *fx = reinterpret_cast<CP64x2 *>(f);
      auto *gx = reinterpret_cast<const CP64x2 *>(g);
      for (u32 i = 0; i != n / 2; ++i)
        fx[i] *= gx[i];
    }
  }
  void div2n(CP64 *f, u32 n) {
    f64 ivn = f64(1) / n;
    if (n <= 16) {
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
