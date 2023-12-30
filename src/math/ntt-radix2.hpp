#ifndef ALGO_H_MATH_NTT_RADIX2
#define ALGO_H_MATH_NTT_RADIX2

#include "../base.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT, u32 G>
struct NttR2 {
  inline static std::array<ModT, 64> rt, irt, rate2, irate2;
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
  }
  static void ntt(void *p, u32 n) {
    auto *f = reinterpret_cast<ModT *>(p);
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
  static void intt(void *p, u32 n) {
    auto *f = reinterpret_cast<ModT *>(p);
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
