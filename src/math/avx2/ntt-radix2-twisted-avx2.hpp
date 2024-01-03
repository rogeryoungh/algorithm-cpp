#ifndef ALGO_H_MATH_AVX2_NTT_RADIX2_TWISTED
#define ALGO_H_MATH_AVX2_NTT_RADIX2_TWISTED

#include "../../base.hpp"
#include "../../modular/avx2/mont32x8.hpp"
#include "../../other/align-alloc.hpp"
#include <algorithm>

ALGO_BEGIN_NAMESPACE

template <class ModT, u32 G>
struct NTT32Radix2TwistedAVX2 {
  using M32x8 = struct M32x8<ModT>;
  // inline static std::array<ModT, 64> rt, irt, rate2, irate2;
  inline static AVec<ModT> rt;
  inline static M32x8 rt2, rt4;
  static void set_mod() {
    rt = {ModT{1}, ModT{1}};
    prepare_root(16);
    rt2 = rotate<4>(rt[3]);
    rt4 = rotate<8>(rt[5]);
  }
  static void prepare_root(u32 m) {
    u32 n = rt.size();
    if (n >= m)
      return;
    rt.resize(m);
    for (; n != m; n *= 2) {
      if (n < 32) {
        ModT w = ModT(G).pow((ModT::MOD - 1) / n / 2);
        for (u32 i = n; i != n * 2; i += 2) {
          rt[i] = rt[i / 2], rt[i + 1] = w * rt[i];
        }
      } else {
        ModT w = ModT(G).pow((ModT::MOD - 1) / n / 2);
        for (u32 i = n; i != n * 2; i += 2) {
          rt[i] = rt[i / 2], rt[i + 1] = w * rt[i];
        }
      }
    }
  }
  static void ntt_butterfly(M32x8 *f, u32 l, M32x8 *w) {
    for (u32 j = 0; j != l; ++j) {
      M32x8 x = f[j], y = f[j + l];
      f[j] = x + y;
      f[j + l] = (x - y) * w[j];
    }
  }
  static void intt_butterfly(M32x8 *f, u32 l, M32x8 *w) {
    for (u32 j = 0; j != l; ++j) {
      M32x8 x = f[j], y = f[j + l] * w[j];
      f[j] = x + y;
      f[j + l] = x - y;
    }
  }
  static void ntt_small(ModT *f, u32 n) {
    for (u32 l = n / 2; l != 0; l /= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        for (u32 j = 0; j != l; ++j) {
          ModT x = f[i + j], y = f[i + j + l];
          f[i + j] = x + y;
          f[i + j + l] = (x - y) * rt[l + j];
        }
      }
    }
  }
  static void intt_small(ModT *f, u32 n) {
    for (u32 l = 1; l != n; l *= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        for (u32 j = 0; j != l; ++j) {
          ModT x = f[i + j], y = f[i + j + l] * rt[l + j];
          f[i + j] = x + y;
          f[i + j + l] = x - y;
        }
      }
    }
  }
  static void ntt_base(M32x8 *f, u32 n, M32x8 *rtx) {
    for (u32 l = n / 2; l != 0; l /= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        ntt_butterfly(f + i, l, rtx + l);
      }
    }
  }
  static void intt_base(M32x8 *f, u32 n, M32x8 *rtx) {
    for (u32 l = 1; l != n; l *= 2) {
      for (u32 i = 0; i != n; i += l * 2) {
        intt_butterfly(f + i, l, rtx + l);
      }
    }
  }
  static void ntt_rec(M32x8 *f, u32 n, M32x8 *rtx) {
    constexpr u32 N = 1 << 6;
    if (n <= N) {
      ntt_base(f, n, rtx);
    } else {
      u32 l = n / 2;
      ntt_butterfly(f, l, rtx + l);
      ntt_rec(f + 0, l, rtx);
      ntt_rec(f + l, l, rtx);
    }
  }
  static void intt_rec(M32x8 *f, u32 n, M32x8 *rtx) {
    constexpr u32 N = 1 << 6;
    if (n <= N) {
      intt_base(f, n, rtx);
    } else {
      u32 l = n / 2;
      intt_rec(f + 0, l, rtx);
      intt_rec(f + l, l, rtx);
      intt_butterfly(f, l, rtx + l);
    }
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
  static void ntt(ModT *p, u32 n) {
    prepare_root(n);
    if (n <= 32)
      return ntt_small(p, n);
    u32 m = n / 8;
    auto *f = reinterpret_cast<M32x8 *>(p);
    auto *rtx = reinterpret_cast<M32x8 *>(rt.data());
    ntt_rec(f, m, rtx);
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = f[i], a, b;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      fi = (a + b) * rt4;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * rt2;
      a = fi.template neg<0xaa>(), b = u32x8_shuffle<0xb1>(fi);
      f[i] = (a + b);
    }
  }
  static void intt(ModT *p, u32 n) {
    prepare_root(n);
    if (n <= 32) {
      intt_small(p, n);
      return std::reverse(p + 1, p + n);
    }
    u32 m = n / 8;
    auto *f = reinterpret_cast<M32x8 *>(p);
    auto *rtx = reinterpret_cast<M32x8 *>(rt.data());
    for (u32 i = 0; i != m; ++i) {
      M32x8 fi = f[i], a, b;
      a = fi.template neg<0xaa>(), b = u32x8_shuffle<0xb1>(fi);
      fi = (a + b) * rt2;
      a = fi.template neg<0xcc>(), b = u32x8_shuffle<0x4e>(fi);
      fi = (a + b) * rt4;
      a = fi.template neg<0xf0>(), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      f[i] = (a + b);
    }
    intt_rec(f, m, rtx);
    std::reverse(p + 1, p + n);
  }
  static void dot(ModT *f, const ModT *g, u32 n) {
    if (n <= 16) {
      for (u32 i = 0; i != n; ++i)
        f[i] *= g[i];
    } else {
      auto *f8 = reinterpret_cast<M32x8 *>(f);
      auto *g8 = reinterpret_cast<const M32x8 *>(g);
      for (u32 i = 0; i != n / 8; ++i)
        f8[i] *= g8[i];
    }
  }
  static void dot2(ModT *f, u32 n) {
    ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
    if (n <= 16) {
      for (u32 i = 0; i != n; ++i)
        f[i] *= ivn;
    } else {
      auto *f8 = reinterpret_cast<M32x8 *>(f);
      M32x8 ivn8 = M32x8::from(ivn);
      for (u32 i = 0; i != n / 8; ++i)
        f8[i] *= ivn8;
    }
  }
};

ALGO_END_NAMESPACE

#endif
