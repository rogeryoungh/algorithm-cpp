#ifndef ALGO_H_MATH_NTT_RADIX4
#define ALGO_H_MATH_NTT_RADIX4

#include "../base.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT, u32 G>
struct NttR4 {
  inline static std::array<ModT, 64> rt, irt, rate2, irate2, rate3, irate3;
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
      rate3[i] = prod * rt[i + 3];
      irate3[i] = iprod * irt[i + 3];
      prod *= irt[i + 3];
      iprod *= rt[i + 3];
    }
  }
  static void ntt(void *p, u32 n) {
    auto *f = reinterpret_cast<ModT *>(p);
    u32 l = n / 2, n4b = std::countr_zero(n) & 1;
    if (n4b) {
      for (u32 j = 0; j != l; ++j) {
        ModT x = f[j], y = f[j + l];
        f[j] = x + y, f[j + l] = x - y;
      }
      l /= 2;
    }
    for (l /= 2; l >= 1; l /= 4) {
      ModT r = 1, img = rt[2];
      for (u32 i = 0, k = 0; i != n; i += l * 4, ++k) {
        ModT r2 = r * r, r3 = r2 * r;
        for (u32 j = 0; j != l; ++j) {
          ModT x0 = f[i + j + 0 * l];
          ModT x1 = f[i + j + 1 * l] * r;
          ModT x2 = f[i + j + 2 * l] * r2;
          ModT x3 = f[i + j + 3 * l] * r3;
          ModT x1x3 = (x1 - x3) * img;
          ModT x02 = x0 + x2, x0_2 = x0 - x2;
          ModT x13 = x1 + x3;
          f[i + j + 0 * l] = x02 + x13;
          f[i + j + 1 * l] = x02 - x13;
          f[i + j + 2 * l] = x0_2 + x1x3;
          f[i + j + 3 * l] = x0_2 - x1x3;
        }
        r *= rate3[std::countr_one(k)];
      }
    }
  }
  static void intt(void *p, u32 n) {
    auto *f = reinterpret_cast<ModT *>(p);
    u32 l = 1, n4b = std::countr_zero(n) & 1;
    for (; l != (n4b ? n / 2 : n); l *= 4) {
      ModT r = 1, img = irt[2];
      for (u32 i = 0, k = 0; i != n; i += l * 4, ++k) {
        ModT r2 = r * r, r3 = r2 * r;
        for (u32 j = 0; j != l; ++j) {
          ModT x0 = f[i + j + 0 * l];
          ModT x1 = f[i + j + 1 * l];
          ModT x2 = f[i + j + 2 * l];
          ModT x3 = f[i + j + 3 * l];
          ModT x2x3 = (x2 - x3) * img;
          ModT x01 = x0 + x1, x0_1 = x0 - x1;
          ModT x23 = x2 + x3;
          f[i + j + 0 * l] = x01 + x23;
          f[i + j + 1 * l] = (x0_1 + x2x3) * r;
          f[i + j + 2 * l] = (x01 - x23) * r2;
          f[i + j + 3 * l] = (x0_1 - x2x3) * r3;
        }
        r *= irate3[std::countr_one(k)];
      }
    }
    if (n4b) {
      for (u32 j = 0; j != l; ++j) {
        ModT x = f[j], y = f[j + l];
        f[j] = x + y, f[j + l] = x - y;
      }
    }
  }
  static void dot(void *p1, const void *p2, u32 n) {
    auto *f = reinterpret_cast<ModT *>(p1);
    auto *g = reinterpret_cast<const ModT *>(p2);
    for (u32 i = 0; i != n; ++i)
      f[i] *= g[i];
  }
  static void dot2(void *p, u32 n) {
    auto *f = reinterpret_cast<ModT *>(p);
    ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
    for (u32 i = 0; i != n; ++i)
      f[i] *= ivn;
  }
};

ALGO_END_NAMESPACE

#endif
