#ifndef ALGO_H_MATH_AVX2_FFT_RADIX2
#define ALGO_H_MATH_AVX2_FFT_RADIX2

#include "./fft-base-avx2.hpp"
#include "../../other/align-alloc.hpp"
#include <numbers>

ALGO_BEGIN_NAMESPACE

struct FFT64Radix2AVX2 : FFT64BaseAVX2 {
  AVec<CP64> rt;
  FFT64Radix2AVX2() {
    rt = {CP64{1}};
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    for (; n != m; n *= 2) {
      CP64 w = CP64::polar(std::numbers::pi / n / 2);
      if (n < 32) {
        for (u32 i = n; i != n * 2; ++i) {
          rt[i] = rt[i - n] * w;
        }
      } else {
        auto *rtx = reinterpret_cast<CP64x2 *>(rt.data());
        CP64x2 wx = CP64x2::from(w);
        for (u32 i = n / 2; i != n; ++i) {
          rtx[i] = rtx[i - n / 2] * wx;
        }
      }
    }
  }
  void fft_small(CP64 *f, u32 n) {
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
  void fft(CP64 *f, u32 n) {
    prepare_root(n);
    if (n <= 8)
      return fft_small(f, n);
    u32 m = n / 2;
    auto *fx = reinterpret_cast<CP64x2 *>(f);
    auto *rtx = reinterpret_cast<CP64x2 *>(rt.data());
    for (u32 l = m / 2; l != 0; l /= 2) {
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        CP64x2 r = CP64x2::from(rt[k]);
        for (u32 j = 0; j != l; ++j) {
          CP64x2 x = fx[i + j], y = fx[i + j + l] * r;
          fx[i + j] = x + y;
          fx[i + j + l] = x - y;
        }
      }
    }
    for (u32 i = 0, k = 0; i != m; i += 2, ++k) {
      CP64x2 abcd = fx[i], efgh = fx[i + 1];
      CP64x2 abef = _mm256_permute2f128_pd(abcd, efgh, 0x20);
      CP64x2 cdgh = _mm256_permute2f128_pd(abcd, efgh, 0x31);
      cdgh *= rtx[k];
      CP64x2 add = abef + cdgh;
      CP64x2 sub = abef - cdgh;
      fx[i] = _mm256_permute2f128_pd(add, sub, 0x20);
      fx[i + 1] = _mm256_permute2f128_pd(add, sub, 0x31);
    }
  }
  void ifft_small(CP64 *f, u32 n) {
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
  void ifft(CP64 *f, u32 n) {
    prepare_root(n);
    if (n <= 8)
      return ifft_small(f, n);
    u32 m = n / 2;
    auto *fx = reinterpret_cast<CP64x2 *>(f);
    auto *rtx = reinterpret_cast<CP64x2 *>(rt.data());
    for (u32 i = 0, k = 0; i != m; i += 2, ++k) {
      CP64x2 abcd = fx[i], efgh = fx[i + 1];
      CP64x2 abef = _mm256_permute2f128_pd(abcd, efgh, 0x20);
      CP64x2 cdgh = _mm256_permute2f128_pd(abcd, efgh, 0x31);
      CP64x2 add = abef + cdgh;
      CP64x2 sub = abef - cdgh;
      sub = CP64x2::cmul(sub, rtx[k]);
      fx[i] = _mm256_permute2f128_pd(add, sub, 0x20);
      fx[i + 1] = _mm256_permute2f128_pd(add, sub, 0x31);
    }
    for (u32 l = 1; l != m; l *= 2) {
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        CP64x2 r = CP64x2::from(rt[k]);
        for (u32 j = 0; j != l; ++j) {
          CP64x2 x = fx[i + j], y = fx[i + j + l];
          fx[i + j] = x + y;
          fx[i + j + l] = CP64x2::cmul(x - y, r);
        }
      }
    }
  }
};

ALGO_END_NAMESPACE

#endif
