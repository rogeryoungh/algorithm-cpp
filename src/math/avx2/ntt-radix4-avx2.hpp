#ifndef ALGO_H_MATH_AVX2_NTT_RADIX2
#define ALGO_H_MATH_AVX2_NTT_RADIX2

#include "../../base.hpp"
#include "../../modular/avx2/mont32x8.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT, u32 G>
struct Ntt32R4Avx2 {
  using M32x8 = struct M32x8<ModT>;
  inline static std::array<ModT, 32> rt, irt, rate2, irate2;
  inline static M32x8 rate3ix8[32], irate3ix8[32], rate4ix8[32], irate4ix8[32];
  inline static M32x8 rt2, irt2, rt4, irt4;
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
    for (u32 i = 0; i != rank2 - 2; ++i) {
      alignas(i256) std::array<ModT, 8> a, ia;
      ModT rate3 = prod * rt[i + 3];
      ModT irate3 = iprod * irt[i + 3];
      for (u32 i = 0; i != 4; ++i) {
        a[i + 4] = a[i] = (i == 0 ? 1 : a[i - 1] * rate3);
        ia[i + 4] = ia[i] = (i == 0 ? 1 : ia[i - 1] * irate3);
      }
      rate3ix8[i] = i256_load(&a);
      irate3ix8[i] = i256_load(&ia);
      prod *= irt[i + 3];
      iprod *= rt[i + 3];
    }

    prod = 1, iprod = 1;
    for (u32 i = 0; i != rank2 - 3; ++i) {
      rate4ix8[i] = rotate(prod * rt[i + 4]);
      irate4ix8[i] = rotate(iprod * irt[i + 4]);
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
    auto *f = reinterpret_cast<M32x8 *>(p);
    u32 m = n / 8, l = m / 2, n4b = std::countr_zero(m) & 1;
    if (n4b) {
      for (u32 j = 0; j != l; ++j) {
        M32x8 x = f[j], y = f[j + l];
        f[j] = x + y, f[j + l] = x - y;
      }
      l /= 2;
    }
    for (l /= 2; l != 0; l /= 4) {
      M32x8 rx = M32x8::from(ModT(1));
      M32x8 rt2x = M32x8::from(rt[2]);
      for (u32 i = 0, k = 0; i != m; i += l * 4, ++k) {
        M32x8 rx8 = u32x8_shuffle<0x55>(rx);
        M32x8 r2x8 = u32x8_shuffle<0xaa>(rx);
        M32x8 r3x8 = u32x8_shuffle<0xff>(rx);
        for (u32 j = 0; j != l; ++j) {
          M32x8 x0 = f[i + j + 0 * l];
          M32x8 x1 = f[i + j + 1 * l];
          M32x8 x2 = f[i + j + 2 * l];
          M32x8 x3 = f[i + j + 3 * l];
          x1 *= rx8, x2 *= r2x8, x3 *= r3x8;
          M32x8 x1x3 = (x1 - x3) * rt2x;
          M32x8 x02 = x0 + x2, x0_2 = x0 - x2;
          M32x8 x13 = x1 + x3;
          f[i + j + 0 * l] = x02 + x13;
          f[i + j + 1 * l] = x02 - x13;
          f[i + j + 2 * l] = x0_2 + x1x3;
          f[i + j + 3 * l] = x0_2 - x1x3;
        }
        rx *= rate3ix8[std::countr_one(k)];
      }
    }
    M32x8 rti = M32x8::from(ModT(1));
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = f[i] * rti, a, b;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      fi = (a + b) * rt4;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * rt2;
      a = fi.template neg<0xaa>(), b = u32x8_shuffle<0xb1>(fi);
      f[i] = (a + b);
      rti *= rate4ix8[std::countr_one(i)];
    }
  }
  static void intt_small(ModT *f, u32 n) {
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
  static void intt(void *p, u32 n) {
    if (n <= 32)
      return intt_small(reinterpret_cast<ModT *>(p), n);
    auto *f = reinterpret_cast<M32x8 *>(p);
    u32 m = n / 8, l = 1, n4b = std::countr_zero(m) & 1;
    M32x8 rti = M32x8::from(ModT(1));
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = f[i], a, b;
      a = fi.template neg<0xaa>(), b = u32x8_shuffle<0xb1>(fi);
      fi = (a + b) * irt2;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * irt4;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      f[i] = (a + b) * rti;
      rti *= irate4ix8[std::countr_one(i)];
    }
    for (; l != (n4b ? m / 2 : m); l *= 4) {
      M32x8 rx = M32x8::from(ModT(1));
      M32x8 irt2x = M32x8::from(irt[2]);
      for (u32 i = 0, k = 0; i != m; i += l * 4, ++k) {
        M32x8 rx8 = u32x8_shuffle<0x55>(rx);
        M32x8 r2x8 = u32x8_shuffle<0xaa>(rx);
        M32x8 r3x8 = u32x8_shuffle<0xff>(rx);
        for (u32 j = 0; j != l; ++j) {
          M32x8 x0 = f[i + j + 0 * l];
          M32x8 x1 = f[i + j + 1 * l];
          M32x8 x2 = f[i + j + 2 * l];
          M32x8 x3 = f[i + j + 3 * l];
          M32x8 x2x3 = (x2 - x3) * irt2x;
          M32x8 x01 = x0 + x1, x0_1 = x0 - x1;
          M32x8 x23 = x2 + x3;
          f[i + j + 0 * l] = x01 + x23;
          f[i + j + 1 * l] = (x0_1 + x2x3) * rx8;
          f[i + j + 2 * l] = (x01 - x23) * r2x8;
          f[i + j + 3 * l] = (x0_1 - x2x3) * r3x8;
        }
        rx *= irate3ix8[std::countr_one(k)];
      }
    }
    if (n4b) {
      for (u32 j = 0; j != l; ++j) {
        M32x8 x = f[j], y = f[j + l];
        f[j] = x + y, f[j + l] = x - y;
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
      auto *f8 = reinterpret_cast<M32x8 *>(p1);
      auto *g8 = reinterpret_cast<const M32x8 *>(p2);
      for (u32 i = 0; i != n / 8; ++i)
        f8[i] *= g8[i];
    }
  }
  static void dot2(void *p, u32 n) {
    ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
    if (n <= 16) {
      auto *f = reinterpret_cast<ModT *>(p);
      for (u32 i = 0; i != n; ++i)
        f[i] *= ivn;
    } else {
      auto *f8 = reinterpret_cast<M32x8 *>(p);
      M32x8 ivn8 = M32x8::from(ivn);
      for (u32 i = 0; i != n / 8; ++i)
        f8[i] *= ivn8;
    }
  }
};

ALGO_END_NAMESPACE

#endif
