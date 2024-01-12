#ifndef ALGO_H_MATH_NTT_RADIX2_TWISTED
#define ALGO_H_MATH_NTT_RADIX2_TWISTED

#include "../base.hpp"
#include "./ntt-base.hpp"
#include <algorithm>
#include <vector>

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTTRadix2Twisted  : NTTBase<ModT> {
  u32 G;
  std::vector<ModT> rt;
  NTTRadix2Twisted(u32 g) : G(g) {
    rt = {ModT{1}, ModT{1}};
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    for (; n != m; n *= 2) {
      ModT w = ModT(G).pow((ModT::MOD - 1) / n / 2);
      for (u32 i = n; i != n * 2; i += 2) {
        rt[i] = rt[i / 2], rt[i + 1] = w * rt[i];
      }
    }
  }
  void ntt_butterfly(ModT *f, u32 l, ModT *w) {
    for (u32 j = 0; j != l; ++j) {
      ModT x = f[j], y = f[j + l];
      f[j] = x + y;
      f[j + l] = (x - y) * w[j];
    }
  }
  void intt_butterfly(ModT *f, u32 l, ModT *w) {
    for (u32 j = 0; j != l; ++j) {
      ModT x = f[j], y = f[j + l] * w[j];
      f[j] = x + y;
      f[j + l] = x - y;
    }
  }
  void ntt_base(ModT *f, u32 n) {
    for (u32 l = n / 2; l != 0; l /= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        ntt_butterfly(f + i, l, rt.data() + l);
      }
    }
  }
  void intt_base(ModT *f, u32 n) {
    for (u32 l = 1; l != n; l *= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        intt_butterfly(f + i, l, rt.data() + l);
      }
    }
  }
  void ntt_rec(ModT *f, u32 n) {
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
  void intt_rec(ModT *f, u32 n) {
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
  void ntt(ModT *f, u32 n) {
    prepare_root(n);
    ntt_base(f, n);
  }
  void intt(ModT *f, u32 n) {
    prepare_root(n);
    intt_base(f, n);
    std::reverse(f + 1, f + n);
  }
};

ALGO_END_NAMESPACE

#endif
