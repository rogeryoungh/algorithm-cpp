#ifndef ALGO_H_MATH_AVX512F_NTT_RADIX2
#define ALGO_H_MATH_AVX512F_NTT_RADIX2

#include "./ntt-base-avx512f.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTT32Radix2AVX512F : NTT32BaseAVX512F<ModT> {
  using M32x16 = struct M32x16<ModT>;
  std::array<ModT, 32> rt, irt, rate2, irate2;
  M32x16 rate5ix[32], irate5ix[32];
  M32x16 rt2x, irt2x, rt4x, irt4x, rt8x, irt8x;
  NTT32Radix2AVX512F(u32 G) {
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
    for (u32 i = 0; i != rank2 - 4; ++i) {
      rate5ix[i] = rotate(prod * rt[i + 5]);
      irate5ix[i] = rotate(iprod * irt[i + 5]);
      prod *= irt[i + 5];
      iprod *= rt[i + 5];
    }
    rt2x = rotate<4>(rt[2]);
    irt2x = rotate<4>(irt[2]);
    rt4x = rotate<8>(rt[3]);
    irt4x = rotate<8>(irt[3]);
    rt8x = rotate<16>(rt[4]);
    irt8x = rotate<16>(irt[4]);
  }

  template <u32 k = 0>
  static M32x16 rotate(ModT t) {
    alignas(i512) std::array<ModT, 16> a;
    if constexpr (k == 0) {
      for (u32 i = 0; i != 16; ++i)
        a[i] = i == 0 ? 1 : a[i - 1] * t;
    } else { // half
      for (u32 i = 0; i != 16; i += k)
        for (u32 j = 0; j != k; ++j)
          a[i + j] = (j <= k / 2) ? 1 : a[i + j - 1] * t;
    }
    return i512_load(&a);
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
    u32 m = n / 16;
    auto *fx = reinterpret_cast<M32x16 *>(f);
    for (u32 l = m / 2; l != 0; l /= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        M32x16 rx = M32x16::from(r);
        for (u32 j = 0; j != l; ++j) {
          M32x16 x = fx[i + j], y = fx[i + j + l] * rx;
          fx[i + j] = x + y;
          fx[i + j + l] = x - y;
        }
        r *= rate2[std::countr_one(k)];
      }
    }
    M32x16 rtix = M32x16::from(ModT(1));
    i512 id1x = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    for (u32 i = 0; i != m; ++i) {
      M32x16 fi = fx[i] * rtix, a, b;
      a = fi.template neg<0xff00>(), b = _mm512_permutexvar_epi64(id1x, fi);
      fi = (a + b) * rt8x;
      a = fi.template neg<0xf0f0>(), b = _mm512_permutex_epi64(fi, 0x4e);
      fi = (a + b) * rt4x;
      a = fi.template neg<0xcccc>(), b = u32x16_shuffle<_MM_PERM_BADC>(fi);
      fi = (a + b) * rt2x;
      a = fi.template neg<0xaaaa>(), b = u32x16_shuffle<_MM_PERM_CDAB>(fi);
      fx[i] = (a + b);
      rtix *= rate5ix[std::countr_one(i)];
    }
  }
  void intt_small(ModT *f, u32 n) {
    for (u32 l = 1; l < n; l *= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i < n; i += l * 2, ++k) {
        for (u32 j = 0; j < l; ++j) {
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
    auto *fx = reinterpret_cast<M32x16 *>(f);
    u32 m = n / 16;
    M32x16 rtix = M32x16::from(ModT(1));
    i512 id1x = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    for (u32 i = 0; i != m; ++i) {
      M32x16 fi = fx[i], a, b;
      a = fi.template neg<0xaaaa>(), b = u32x16_shuffle<_MM_PERM_CDAB>(fi);
      fi = (a + b) * irt2x;
      a = fi.template neg<0xcccc>(), b = u32x16_shuffle<_MM_PERM_BADC>(fi);
      fi = (a + b) * irt4x;
      a = fi.template neg<0xf0f0>(), b = _mm512_permutex_epi64(fi, 0x4e);
      fi = (a + b) * irt8x;
      a = fi.template neg<0xff00>(), b = _mm512_permutexvar_epi64(id1x, fi);
      fx[i] = (a + b) * rtix;
      rtix *= irate5ix[std::countr_one(i)];
    }
    for (u32 l = 1; l != m; l *= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        M32x16 rx = M32x16::from(r);
        for (u32 j = 0; j != l; ++j) {
          M32x16 x = fx[i + j], y = fx[i + j + l];
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
