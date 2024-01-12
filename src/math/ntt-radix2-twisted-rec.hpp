#ifndef ALGO_H_MATH_NTT_RADIX2_TWISTED_REC
#define ALGO_H_MATH_NTT_RADIX2_TWISTED_REC

#include "../base.hpp"
#include "./ntt-base.hpp"
#include <vector>
#include <iostream>

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTTRadix2Twisted : NTTBase<ModT> {
  u32 G;
  std::vector<ModT> rt, irt;
  NTTRadix2Twisted(u32 g) : G(g) {
    irt = rt = {ModT{1}, ModT{1}};
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    irt.resize(m);
    for (; n != m; n *= 2) {
      ModT w = ModT(G).pow((ModT::MOD - 1) / n / 2);
      ModT iw = w.inv();
      for (u32 i = n; i != n * 2; i += 2) {
        rt[i] = rt[i / 2], rt[i + 1] = w * rt[i];
        irt[i] = irt[i / 2], irt[i + 1] = iw * irt[i];
      }
    }
  }
  template <u32 n>
  void ntt_rec_t(ModT *f) {
    if constexpr (n <= 1) {
      return;
    } else if constexpr (n == 2) {
      ModT f0 = f[0];
      ModT f1 = f[1];
      f[0] = f0 + f1;
      f[1] = f0 - f1;
    } else {
      for (u32 k = 0; k != n / 2; ++k) {
        ModT f0 = f[k];
        ModT f1 = f[k + n / 2];
        f[k] = f0 + f1;
        f[k + n / 2] = (f0 - f1) * irt[n / 2 + k];
      }
      ntt_rec_t<n / 2>(f);
      ntt_rec_t<n / 2>(f + n / 2);
    }
  }
  template <u32 n>
  void intt_rec_t(ModT *f) {
    if constexpr (n <= 1) {
      return;
    } else if constexpr (n == 2) {
      ModT f0 = f[0];
      ModT f1 = f[1];
      f[0] = f0 + f1;
      f[1] = f0 - f1;
    } else {
      intt_rec_t<n / 2>(f);
      intt_rec_t<n / 2>(f + n / 2);
      for (u32 k = 0; k != n / 2; ++k) {
        ModT f0 = f[k];
        ModT f1 = f[k + n / 2] * rt[n / 2 + k];
        f[k] = f0 + f1;
        f[k + n / 2] = f0 - f1;
      }
    }
  }
  template <u32 n>
  void ntt_rec(ModT *f, u32 len) {
    if constexpr (n > 0) {
      if (n == len) {
        ntt_rec_t<n>(f);
      } else {
        ntt_rec<n * 2>(f, len);
      }
    }
  }
  template <u32 n>
  void intt_rec(ModT *f, u32 len) {
    if constexpr (n > 0) {
      if (n == len) {
        intt_rec_t<n>(f);
      } else {
        intt_rec<n * 2>(f, len);
      }
    }
  }
  void ntt(ModT *f, u32 n) {
    prepare_root(n);
    ntt_rec<1>(f, n);
  }
  void intt(ModT *f, u32 n) {
    prepare_root(n);
    intt_rec<1>(f, n);
  }
};

ALGO_END_NAMESPACE

#endif
