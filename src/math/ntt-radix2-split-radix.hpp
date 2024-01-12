#ifndef ALGO_H_MATH_NTT_RADIX2_SPLIT_RADIX
#define ALGO_H_MATH_NTT_RADIX2_SPLIT_RADIX

#include "../base.hpp"
#include "./ntt-base.hpp"
#include <vector>

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTTRadix2Split : NTTBase<ModT> {
  u32 G;
  ModT J;
  std::vector<ModT> rt, irt;
  NTTRadix2Split(u32 g) : G(g) {
    irt = rt = {ModT{1}, ModT{1}};
    J = ModT(G).pow((ModT::MOD - 1) / 4);
  }
  void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m), irt.resize(m);
    for (; n != m; n *= 2) {
      ModT w = ModT(G).pow((ModT::MOD - 1) / n / 4);
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
    } else if constexpr (n == 4) {
      for (u32 k = 0; k != n / 4; ++k) {
        ModT f0 = f[k];
        ModT f2 = f[k + n / 2];
        ModT f1 = f[k + n / 4];
        ModT f3 = f[k + 3 * n / 4];
        ModT f02s = f0 - f2;
        ModT f13s = (f1 - f3) * J;
        f[k] = f0 + f2;
        f[k + n / 2] = (f02s - f13s);
        f[k + n / 4] = f1 + f3;
        f[k + 3 * n / 4] = (f02s + f13s);
      }
      ntt_rec_t<n / 2>(f);
    } else {
      for (u32 k = 0; k != n / 4; ++k) {
        ModT f0 = f[k];
        ModT f2 = f[k + n / 2];
        ModT f1 = f[k + n / 4];
        ModT f3 = f[k + 3 * n / 4];
        ModT f02s = f0 - f2;
        ModT f13s = (f1 - f3) * J;
        f[k] = f0 + f2;
        f[k + n / 2] = (f02s - f13s) * irt[n / 4 + k];
        f[k + n / 4] = f1 + f3;
        f[k + 3 * n / 4] = (f02s + f13s) * rt[n / 4 + k];
      }
      ntt_rec_t<n / 2>(f);
      ntt_rec_t<n / 4>(f + n / 2);
      ntt_rec_t<n / 4>(f + 3 * n / 4);
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
    } else if constexpr (n == 4) {
      intt_rec_t<n / 2>(f);
      for (u32 k = 0; k != n / 4; ++k) {
        ModT f0 = f[k];
        ModT f2 = f[k + n / 2];
        ModT f1 = f[k + n / 4];
        ModT f3 = f[k + 3 * n / 4];
        ModT f23a = f2 + f3;
        ModT f23s = (f2 - f3) * J;
        f[k] = f0 + f23a;
        f[k + n / 2] = f0 - f23a;
        f[k + n / 4] = f1 + f23s;
        f[k + 3 * n / 4] = f1 - f23s;
      }
    } else {
      intt_rec_t<n / 2>(f);
      intt_rec_t<n / 4>(f + n / 2);
      intt_rec_t<n / 4>(f + 3 * n / 4);
      for (u32 k = 0; k != n / 4; ++k) {
        ModT f0 = f[k];
        ModT f2 = f[k + n / 2] * rt[n / 4 + k];
        ModT f1 = f[k + n / 4];
        ModT f3 = f[k + 3 * n / 4] * irt[n / 4 + k];
        ModT f23a = f2 + f3;
        ModT f23s = (f2 - f3) * J;
        f[k] = f0 + f23a;
        f[k + n / 2] = f0 - f23a;
        f[k + n / 4] = f1 + f23s;
        f[k + 3 * n / 4] = f1 - f23s;
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
    prepare_root(n / 2);
    ntt_rec<1>(f, n);
  }
  void intt(ModT *f, u32 n) {
    prepare_root(n / 2);
    intt_rec<1>(f, n);
  }
};

ALGO_END_NAMESPACE

#endif
