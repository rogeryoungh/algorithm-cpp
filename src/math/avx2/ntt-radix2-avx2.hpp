#ifndef ALGO_H_MATH_AVX2_NTT_RADIX2
#define ALGO_H_MATH_AVX2_NTT_RADIX2

#include "./ntt-base-avx2.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTT32Radix2AVX2 : NTT32BaseAVX2<ModT> {
  using M32x8 = struct M32x8<ModT>;
  std::array<ModT, 32> rt, irt, rate2, irate2;
  M32x8 rate4ix[32], irate4ix[32];
  M32x8 rt2, irt2, rt4, irt4;
  NTT32Radix2AVX2(u32 G) {
    const u32 rank2 = std::countr_zero(ModT::MOD - 1);
    rt[rank2] = ModT(G).pow((ModT::MOD - 1) >> rank2);
    irt[rank2] = rt[rank2].inv();
    for (u32 i = rank2; i != 1; --i) {
      rt[i - 1] = rt[i] * rt[i];
      irt[i - 1] = irt[i] * irt[i];
    }
    ModT prod = 1, iprod = 1;
    for (u32 i = 0; i != rank2 - 1; ++i) {
      rate2[i] = prod * rt[i + 2];
      irate2[i] = iprod * irt[i + 2];
      prod *= irt[i + 2];
      iprod *= rt[i + 2];
    }

    prod = 1, iprod = 1;
    for (u32 i = 0; i != rank2 - 3; ++i) {
      rate4ix[i] = rotate(prod * rt[i + 4]);
      irate4ix[i] = rotate(iprod * irt[i + 4]);
      prod *= irt[i + 4];
      iprod *= rt[i + 4];
    }
    rt2 = rotate<4>(rt[2]);
    irt2 = rotate<4>(irt[2]);
    rt4 = rotate<8>(rt[3]);
    irt4 = rotate<8>(irt[3]);
  }

  template <u32 k = 0>
  static M32x8 rotate(ModT t) {
    alignas(i256) std::array<ModT, 8> a;
    if constexpr (k == 0) {
      for (u32 i = 0; i != 8; ++i)
        a[i] = i == 0 ? 1 : a[i - 1] * t;
    } else { // half
      for (u32 i = 0; i != 8; i += k)
        for (u32 j = 0; j != k; ++j)
          a[i + j] = (j <= k / 2) ? 1 : a[i + j - 1] * t;
    }
    return i256_load(&a);
  }
  void ntt_small(ModT *f, u32 n) {
    for (u32 l = n / 2; l != 0; l /= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          ModT x = f[i + j], y = f[i + j + l] * r;
          f[i + j] = x + y;
          f[i + j + l] = x - y;
        }
        r *= rate2[std::countr_one(k)];
      }
    }
  }
  void ntt(ModT *f, u32 n) {
    if (n <= 32)
      return ntt_small(f, n);
    u32 m = n / 8;
    auto *fx = reinterpret_cast<M32x8 *>(f);
    for (u32 l = m / 2; l != 0; l /= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        M32x8 rx = M32x8::from(r);
        for (u32 j = 0; j != l; ++j) {
          M32x8 x = fx[i + j], y = fx[i + j + l] * rx;
          fx[i + j] = x + y;
          fx[i + j + l] = x - y;
        }
        r *= rate2[std::countr_one(k)];
      }
    }
    M32x8 rtix = M32x8::from(ModT(1));
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = fx[i] * rtix, a, b;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      fi = (a + b) * rt4;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * rt2;
      a = fi.template neg<0xaa>(), b = u32x8_shuffle<0xb1>(fi);
      fx[i] = (a + b);
      rtix *= rate4ix[std::countr_one(i)];
    }
  }
  void intt_small(ModT *f, u32 n) {
    for (u32 l = 1; l != n; l *= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          ModT x = f[i + j], y = f[i + j + l];
          f[i + j] = x + y;
          f[i + j + l] = (x - y) * r;
        }
        r *= irate2[std::countr_one(k)];
      }
    }
  }
  void intt(ModT *f, u32 n) {
    if (n <= 32)
      return intt_small(f, n);
    u32 m = n / 8;
    auto *fx = reinterpret_cast<M32x8 *>(f);
    M32x8 rti = M32x8::from(ModT(1));
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = fx[i], a, b;
      a = fi.template neg<0xaa>(), b = u32x8_shuffle<0xb1>(fi);
      fi = (a + b) * irt2;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * irt4;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      fx[i] = (a + b) * rti;
      rti *= irate4ix[std::countr_one(i)];
    }
    for (u32 l = 1; l != m; l *= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        M32x8 rx = M32x8::from(r);
        for (u32 j = 0; j != l; ++j) {
          M32x8 x = fx[i + j], y = fx[i + j + l];
          fx[i + j] = x + y;
          fx[i + j + l] = (x - y) * rx;
        }
        r *= irate2[std::countr_one(k)];
      }
    }
  }
};

ALGO_END_NAMESPACE

#endif
