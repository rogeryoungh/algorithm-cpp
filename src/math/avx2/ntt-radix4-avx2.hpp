#ifndef ALGO_H_MATH_AVX2_NTT_RADIX4
#define ALGO_H_MATH_AVX2_NTT_RADIX4

#include "./ntt-base-avx2.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTT32Radix4AVX2 : NTT32BaseAVX2<ModT> {
  using M32x8 = struct M32x8<ModT>;
  std::array<ModT, 32> rt, irt, rate2, irate2;
  M32x8 rate3ix[32], irate3ix[32], rate4ix[32], irate4ix[32];
  M32x8 rt2x, irt2x, rt4x, irt4x;
  NTT32Radix4AVX2(u32 G) {
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
    for (u32 i = 0; i != rank2 - 2; ++i) {
      alignas(i256) std::array<ModT, 8> a, ia;
      ModT rate3 = prod * rt[i + 3];
      ModT irate3 = iprod * irt[i + 3];
      for (u32 i = 0; i != 4; ++i) {
        a[i + 4] = a[i] = (i == 0 ? 1 : a[i - 1] * rate3);
        ia[i + 4] = ia[i] = (i == 0 ? 1 : ia[i - 1] * irate3);
      }
      rate3ix[i] = i256_load(&a);
      irate3ix[i] = i256_load(&ia);
      prod *= irt[i + 3];
      iprod *= rt[i + 3];
    }

    prod = 1, iprod = 1;
    for (u32 i = 0; i != rank2 - 3; ++i) {
      rate4ix[i] = rotate(prod * rt[i + 4]);
      irate4ix[i] = rotate(iprod * irt[i + 4]);
      prod *= irt[i + 4];
      iprod *= rt[i + 4];
    }
    rt2x = rotate<4>(rt[2]);
    irt2x = rotate<4>(irt[2]);
    rt4x = rotate<8>(rt[3]);
    irt4x = rotate<8>(irt[3]);
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
    auto *fx = reinterpret_cast<M32x8 *>(f);
    u32 m = n / 8, l = m / 2, n4b = std::countr_zero(m) & 1;
    if (n4b) {
      for (u32 j = 0; j != l; ++j) {
        M32x8 x = fx[j], y = fx[j + l];
        fx[j] = x + y, fx[j + l] = x - y;
      }
      l /= 2;
    }
    for (l /= 2; l != 0; l /= 4) {
      M32x8 rx = M32x8::from(ModT(1));
      M32x8 z4x = M32x8::from(rt[2]);
      for (u32 i = 0, k = 0; i != m; i += l * 4, ++k) {
        M32x8 r1x = u32x8_shuffle<0x55>(rx);
        M32x8 r2x = u32x8_shuffle<0xaa>(rx);
        M32x8 r3x = u32x8_shuffle<0xff>(rx);
        for (u32 j = 0; j != l; ++j) {
          M32x8 x0 = fx[i + j + 0 * l];
          M32x8 x1 = fx[i + j + 1 * l];
          M32x8 x2 = fx[i + j + 2 * l];
          M32x8 x3 = fx[i + j + 3 * l];
          x1 *= r1x, x2 *= r2x, x3 *= r3x;
          M32x8 x1x3 = (x1 - x3) * z4x;
          M32x8 x02 = x0 + x2, x0_2 = x0 - x2;
          M32x8 x13 = x1 + x3;
          fx[i + j + 0 * l] = x02 + x13;
          fx[i + j + 1 * l] = x02 - x13;
          fx[i + j + 2 * l] = x0_2 + x1x3;
          fx[i + j + 3 * l] = x0_2 - x1x3;
        }
        rx *= rate3ix[std::countr_one(k)];
      }
    }
    M32x8 rtix = M32x8::from(ModT(1));
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = fx[i] * rtix, a, b;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      fi = (a + b) * rt4x;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * rt2x;
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
    auto *fx = reinterpret_cast<M32x8 *>(f);
    u32 m = n / 8, l = 1, n4b = std::countr_zero(m) & 1;
    M32x8 rtix = M32x8::from(ModT(1));
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = fx[i], a, b;
      a = fi.template neg<0xaa>(), b = u32x8_shuffle<0xb1>(fi);
      fi = (a + b) * irt2x;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * irt4x;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      fx[i] = (a + b) * rtix;
      rtix *= irate4ix[std::countr_one(i)];
    }
    for (; l != (n4b ? m / 2 : m); l *= 4) {
      M32x8 rx = M32x8::from(ModT(1));
      M32x8 iz8x = M32x8::from(irt[2]);
      for (u32 i = 0, k = 0; i != m; i += l * 4, ++k) {
        M32x8 r1x8 = u32x8_shuffle<0x55>(rx);
        M32x8 r2x8 = u32x8_shuffle<0xaa>(rx);
        M32x8 r3x8 = u32x8_shuffle<0xff>(rx);
        for (u32 j = 0; j != l; ++j) {
          M32x8 x0 = fx[i + j + 0 * l];
          M32x8 x1 = fx[i + j + 1 * l];
          M32x8 x2 = fx[i + j + 2 * l];
          M32x8 x3 = fx[i + j + 3 * l];
          M32x8 x2x3 = (x2 - x3) * iz8x;
          M32x8 x01 = x0 + x1, x0_1 = x0 - x1;
          M32x8 x23 = x2 + x3;
          fx[i + j + 0 * l] = x01 + x23;
          fx[i + j + 1 * l] = (x0_1 + x2x3) * r1x8;
          fx[i + j + 2 * l] = (x01 - x23) * r2x8;
          fx[i + j + 3 * l] = (x0_1 - x2x3) * r3x8;
        }
        rx *= irate3ix[std::countr_one(k)];
      }
    }
    if (n4b) {
      for (u32 j = 0; j != l; ++j) {
        M32x8 x = fx[j], y = fx[j + l];
        fx[j] = x + y, fx[j + l] = x - y;
      }
    }
  }
};

ALGO_END_NAMESPACE

#endif
