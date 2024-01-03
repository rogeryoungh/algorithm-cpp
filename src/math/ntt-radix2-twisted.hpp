#ifndef ALGO_H_MATH_NTT_RADIX2_TWISTED
#define ALGO_H_MATH_NTT_RADIX2_TWISTED

#include "../base.hpp"
#include <algorithm>
#include <vector>

ALGO_BEGIN_NAMESPACE

template <class ModT, u32 G>
struct NttR2T {
  // inline static std::array<ModT, 64> rt, irt, rate2, irate2;
  inline static std::vector<ModT> rt;
  static void setMod() {
    rt = {ModT{1}, ModT{1}};
  }
  static void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n < m) {
      rt.resize(m);
      for (; n != m; n *= 2) {
        ModT w = ModT(G).pow((ModT::MOD - 1) / n / 2);
        for (u32 i = n; i != n * 2; i += 2) {
          rt[i] = rt[i / 2], rt[i + 1] = w * rt[i];
        }
      }
    }
  }
  static void ntt_butterfly(ModT *f, u32 l, ModT *w) {
    for (u32 j = 0; j != l; ++j) {
      ModT x = f[j], y = f[j + l];
      f[j] = x + y;
      f[j + l] = (x - y) * w[j];
    }
  }
  static void intt_butterfly(ModT *f, u32 l, ModT *w) {
    for (u32 j = 0; j != l; ++j) {
      ModT x = f[j], y = f[j + l] * w[j];
      f[j] = x + y;
      f[j + l] = x - y;
    }
  }
  static void ntt_base(ModT *f, u32 n) {
    for (u32 l = n / 2; l != 0; l /= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        ntt_butterfly(f + i, l, rt.data() + l);
      }
    }
  }
  static void intt_base(ModT *f, u32 n) {
    for (u32 l = 1; l != n; l *= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        intt_butterfly(f + i, l, rt.data() + l);
      }
    }
  }
  static void ntt_rec(ModT *f, u32 n) {
    constexpr u32 N = 1 << 10;
    if (n <= N) {
      ntt_base(f, n);
    } else {
      u32 l = n / 2;
      ntt_butterfly(f, l, rt.data() + l);
      ntt_rec(f + 0, l);
      ntt_rec(f + l, l);
    }
  }
  static void intt_rec(ModT *f, u32 n) {
    constexpr u32 N = 1 << 10;
    if (n <= N) {
      intt_base(f, n);
    } else {
      u32 l = n / 2;
      intt_rec(f + 0, l);
      intt_rec(f + l, l);
      intt_butterfly(f, l, rt.data() + l);
    }
  }
  static void ntt(void *p, u32 n) {
    auto *f = reinterpret_cast<ModT *>(p);
    prepare_root(n);
    ntt_base(f, n);
  }
  static void intt(void *p, u32 n) {
    prepare_root(n);
    auto *f = reinterpret_cast<ModT *>(p);
    intt_base(f, n);
    std::reverse(f + 1, f + n);
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
