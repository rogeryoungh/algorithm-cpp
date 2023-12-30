#ifndef ALGO_H_MATH_AVX512F_NTT_RADIX2
#define ALGO_H_MATH_AVX512F_NTT_RADIX2

#include "../../base.hpp"
#include "../../modular/avx512f/mont32x16.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT, u32 G>
struct Ntt32R2Avx512F {
  using M32x16 = struct M32x16<ModT>;
  inline static std::array<ModT, 32> rt, irt, rate2, irate2;
  inline static M32x16 rate5ix16[32], irate5ix16[32];
  inline static M32x16 rt2, irt2, rt4, irt4, rt8, irt8;
  static void setMod() {
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
      rate5ix16[i] = rotate(prod * rt[i + 5]);
      irate5ix16[i] = rotate(iprod * irt[i + 5]);
      prod *= irt[i + 5];
      iprod *= rt[i + 5];
    }
    rt2 = rotate<4>(rt[2]);
    irt2 = rotate<4>(irt[2]);
    rt4 = rotate<8>(rt[3]);
    irt4 = rotate<8>(irt[3]);
    rt8 = rotate<16>(rt[4]);
    irt8 = rotate<16>(irt[4]);
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
  static void ntt_small(ModT *f, u32 n) {
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
  static void ntt(void *p, u32 n) {
    if (n <= 32)
      return ntt_small(reinterpret_cast<ModT *>(p), n);
    u32 m = n / 16;
    auto *f = reinterpret_cast<M32x16 *>(p);
    for (u32 l = m / 2; l != 0; l /= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        M32x16 rx16 = M32x16::from(r);
        for (u32 j = 0; j != l; ++j) {
          M32x16 x = f[i + j], y = f[i + j + l] * rx16;
          f[i + j] = x + y;
          f[i + j + l] = x - y;
        }
        r *= rate2[std::countr_one(k)];
      }
    }
    M32x16 rti = M32x16::from(ModT(1));
    i512 id1 = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    for (u32 i = 0; i != m; ++i) {
      M32x16 fi = f[i] * rti, a, b;
      a = fi.template neg<0xff00>(), b = _mm512_permutexvar_epi64(id1, fi);
      fi = (a + b) * rt8;
      a = fi.template neg<0xf0f0>(), b = _mm512_permutex_epi64(fi, 0x4e);
      fi = (a + b) * rt4;
      a = fi.template neg<0xcccc>(), b = u32x16_shuffle<_MM_PERM_BADC>(fi);
      fi = (a + b) * rt2;
      a = fi.template neg<0xaaaa>(), b = u32x16_shuffle<_MM_PERM_CDAB>(fi);
      f[i] = (a + b);
      rti *= rate5ix16[std::countr_one(i)];
    }
  }
  static void intt_small(ModT *f, u32 n) {
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
  static void intt(void *p, u32 n) {
    if (n <= 32)
      return intt_small(reinterpret_cast<ModT *>(p), n);
    auto *f = reinterpret_cast<M32x16 *>(p);
    u32 m = n / 16;
    M32x16 rti = M32x16::from(ModT(1));
    i512 id1 = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    for (u32 i = 0; i != m; ++i) {
      M32x16 fi = f[i], a, b;
      a = fi.template neg<0xaaaa>(), b = u32x16_shuffle<_MM_PERM_CDAB>(fi);
      fi = (a + b) * irt2;
      a = fi.template neg<0xcccc>(), b = u32x16_shuffle<_MM_PERM_BADC>(fi);
      fi = (a + b) * irt4;
      a = fi.template neg<0xf0f0>(), b = _mm512_permutex_epi64(fi, 0x4e);
      fi = (a + b) * irt8;
      a = fi.template neg<0xff00>(), b = _mm512_permutexvar_epi64(id1, fi);
      f[i] = (a + b) * rti;
      rti *= irate5ix16[std::countr_one(i)];
    }
    for (u32 l = 1; l != m; l *= 2) {
      ModT r = 1;
      for (u32 i = 0, k = 0; i != m; i += l * 2, ++k) {
        M32x16 rx16 = M32x16::from(r);
        for (u32 j = 0; j != l; ++j) {
          M32x16 x = f[i + j], y = f[i + j + l];
          f[i + j] = x + y;
          f[i + j + l] = (x - y) * rx16;
        }
        r *= irate2[std::countr_one(k)];
      }
    }
  }
  static void dot(void *p1, const void *p2, u32 n) {
    if (n <= 16) {
      auto *f = reinterpret_cast<ModT *>(p1);
      auto *g = reinterpret_cast<const ModT *>(p2);
      for (u32 i = 0; i != n; ++i)
        f[i] *= g[i];
    } else {
      auto *f8 = reinterpret_cast<M32x16 *>(p1);
      auto *g8 = reinterpret_cast<const M32x16 *>(p2);
      for (u32 i = 0; i != n / 16; ++i)
        f8[i] *= g8[i];
    }
  }
  static void dot2(void *p, u32 n) {
    ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
    if (n <= 32) {
      auto *f = reinterpret_cast<ModT *>(p);
      for (u32 i = 0; i != n; ++i)
        f[i] *= ivn;
    } else {
      auto *f16 = reinterpret_cast<M32x16 *>(p);
      M32x16 ivn16 = M32x16::from(ivn);
      for (u32 i = 0; i != n / 16; ++i)
        f16[i] *= ivn16;
    }
  }
};

ALGO_END_NAMESPACE

#endif
