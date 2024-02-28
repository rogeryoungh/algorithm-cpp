#pragma once

#include "../../number/mont32x16.hpp"
#include <array>
#include <bit>

ALGO_BEGIN_NAMESPACE

struct NTT32OriginalRadix2AVX512F {
  std::array<u32, 32> rt, irt, rate2, irate2;
  u32x16 rate5ix[32], irate5ix[32];
  u32x16 _rt2, _irt2, _rt4, _irt4, _rt8, _irt8;
  Mont32x16 _MX;
  NTT32OriginalRadix2AVX512F(Mont32 M, u32 G) : _MX(M) {
    const u32 rank2 = std::countr_zero(M.MOD - 1);
    rt[rank2] = M.qpow(M.trans(G), (M.MOD - 1) >> rank2);
    irt[rank2] = M.inv(rt[rank2]);
    for (u32 i = rank2; i != 0; --i) {
      rt[i - 1] = M.mul(rt[i], rt[i]);
      irt[i - 1] = M.mul(irt[i], irt[i]);
    }
    u32 prod = M.ONE, iprod = M.ONE;
    for (u32 i = 0; i != rank2 - 1; ++i) {
      rate2[i] = M.mul(prod, rt[i + 2]);
      irate2[i] = M.mul(iprod, irt[i + 2]);
      prod = M.mul(prod, irt[i + 2]);
      iprod = M.mul(iprod, rt[i + 2]);
    }
    prod = M.ONE, iprod = M.ONE;
    alignas(i512) u32 arr[16];
    auto rotate = [&M, &arr](u32 x) {
      for (u32 i = 0; i != 16; ++i)
        arr[i] = i == 0 ? M.ONE : M.mul(x, arr[i - 1]);
    };
    for (u32 i = 0; i != rank2 - 4; ++i) {
      rotate(M.mul(prod, rt[i + 5]));
      rate5ix[i] = _MX.loadu(arr);
      rotate(M.mul(iprod, irt[i + 5]));
      irate5ix[i] = _MX.loadu(arr);
      prod = M.mul(prod, irt[i + 5]);
      iprod = M.mul(iprod, rt[i + 5]);
    }
    auto rotatex = [&M, &arr](u32 x, u32 k) {
      for (u32 i = 0; i != 16; i += k)
        for (u32 j = 0; j != k; ++j)
          arr[i + j] = (j <= k / 2) ? M.ONE : M.mul(x, arr[i + j - 1]);
    };
    rotatex(rt[2], 4), _rt2 = _MX.loadu(arr);
    rotatex(irt[2], 4), _irt2 = _MX.loadu(arr);
    rotatex(rt[3], 8), _rt4 = _MX.loadu(arr);
    rotatex(irt[3], 8), _irt4 = _MX.loadu(arr);
    rotatex(rt[4], 16), _rt8 = _MX.loadu(arr);
    rotatex(irt[4], 16), _irt8 = _MX.loadu(arr);
  }
  void ntt_small(u32 *f, usize n) {
    const auto M = _MX.M;
    for (u32 l = n / 2; l >= 1; l /= 2) {
      u32 r = M.ONE;
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        u32 *fx = f + i, *fy = fx + l;
        for (u32 j = 0; j != l; ++j) {
          u32 x = fx[j], y = M.mul(fy[j], r);
          fx[j] = M.add(x, y);
          fy[j] = M.sub(x, y);
        }
        r = M.mul(r, rate2[std::countr_one(k)]);
      }
    }
  }
  void intt_small(u32 *f, usize n) {
    const auto M = _MX.M;
    u32 ivn = M.trans(M.MOD - (M.MOD - 1) / n);
    for (u32 l = 1; l <= n / 2; l *= 2) {
      u32 r = M.ONE;
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        u32 *fx = f + i, *fy = fx + l;
        for (u32 j = 0; j != l; ++j) {
          u32 x = fx[j], y = fy[j];
          fx[j] = M.add(x, y);
          fy[j] = M.mul(M.sub(x, y), r);
        }
        r = M.mul(r, irate2[std::countr_one(k)]);
      }
    }
    for (u32 i = 0; i != n; ++i)
      f[i] = M.mul(f[i], ivn);
  }
  void ntt(u32 *f, usize n) {
    if (n < 16)
      return ntt_small(f, n);
    const auto MX = _MX;
    const auto M = MX.M;
    for (u32 l = n / 2; l >= 16; l /= 2) {
      u32 *f0 = f, *f1 = f + l;
      for (u32 j = 0; j != l; j += 16) {
        u32x16 x = MX.loadu(f0 + j), y = MX.loadu(f1 + j);
        MX.storeu(f0 + j, MX.add(x, y));
        MX.storeu(f1 + j, MX.sub(x, y));
      }
      u32 r = rate2[0];
      for (u32 i = l * 2, k = 1; i != n; i += l * 2, ++k) {
        u32x16 rx = MX.set1(r);
        f0 = f + i, f1 = f0 + l;
        for (u32 j = 0; j != l; j += 16) {
          u32x16 x = MX.loadu(f0 + j), y = MX.mul(rx, MX.loadu(f1 + j));
          MX.storeu(f0 + j, MX.add(x, y));
          MX.storeu(f1 + j, MX.sub(x, y));
        }
        r = M.mul(r, rate2[std::countr_one(k)]);
      }
    }
    u32x16 rtix = MX.ONE, rt2 = _rt2, rt4 = _rt4, rt8 = _rt8;
    i512 id1x = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    for (u32 i = 0; i != n; i += 16) {
      u32x16 fi = MX.mul(rtix, MX.loadu(f + i)), a, b;
      a = MX.neg<0xff00>(fi), b = _mm512_permutexvar_epi64(id1x, fi);
      fi = MX.mul(rt8, MX.add(a, b));
      a = MX.neg<0xf0f0>(fi), b = _mm512_permutex_epi64(fi, 0x4e);
      fi = MX.mul(rt4, MX.add(a, b));
      a = MX.neg<0xcccc>(fi), b = u32x16_shuffle<_MM_PERM_BADC>(fi);
      fi = MX.mul(rt2, MX.add(a, b));
      a = MX.neg<0xaaaa>(fi), b = u32x16_shuffle<_MM_PERM_CDAB>(fi);
      MX.storeu(f + i, MX.add(a, b));
      rtix = MX.mul(rtix, rate5ix[std::countr_one(i / 16)]);
    }
  }
  void intt(u32 *f, usize n) {
    if (n < 16)
      return intt_small(f, n);
    const auto MX = _MX;
    const auto M = MX.M;
    u32x16 rtix = MX.set1(M.trans(M.MOD - (M.MOD - 1) / n));
    u32x16 irt2 = _irt2, irt4 = _irt4, irt8 = _irt8;
    i512 id1x = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    for (u32 i = 0; i != n; i += 16) {
      u32x16 fi = MX.loadu(f + i), a, b;
      a = MX.neg<0xaaaa>(fi), b = u32x16_shuffle<_MM_PERM_CDAB>(fi);
      fi = MX.mul(irt2, MX.add(a, b));
      a = MX.neg<0xcccc>(fi), b = u32x16_shuffle<_MM_PERM_BADC>(fi);
      fi = MX.mul(irt4, MX.add(a, b));
      a = MX.neg<0xf0f0>(fi), b = _mm512_permutex_epi64(fi, 0x4e);
      fi = MX.mul(irt8, MX.add(a, b));
      a = MX.neg<0xff00>(fi), b = _mm512_permutexvar_epi64(id1x, fi);
      MX.storeu(f + i, MX.mul(MX.add(a, b), rtix));
      rtix = MX.mul(rtix, irate5ix[std::countr_one(i / 16)]);
    }
    for (u32 l = 16; l <= n / 2; l *= 2) {
      u32 *f0 = f, *f1 = f + l;
      for (u32 j = 0; j != l; j += 16) {
        u32x16 x = MX.loadu(f0 + j), y = MX.loadu(f1 + j);
        MX.storeu(f0 + j, MX.add(x, y));
        MX.storeu(f1 + j, MX.sub(x, y));
      }
      u32 r = irate2[0];
      for (u32 i = l * 2, k = 1; i != n; i += l * 2, ++k) {
        u32x16 rx = MX.set1(r);
        f0 = f + i, f1 = f0 + l;
        for (u32 j = 0; j != l; j += 16) {
          u32x16 x = MX.loadu(f0 + j), y = MX.loadu(f1 + j);
          MX.storeu(f0 + j, MX.add(x, y));
          MX.storeu(f1 + j, MX.mul(MX.sub(x, y), rx));
        }
        r = M.mul(r, irate2[std::countr_one(k)]);
      }
    }
  }
  void conv(u32 *f, u32 *g, u32 n) {
    if (n < 16) {
      const auto M = _MX.M;
      for (u32 i = 0; i != n; ++i)
        f[i] = M.trans(f[i]), g[i] = M.trans(g[i]);
      ntt(f, n), ntt(g, n);
      for (u32 i = 0; i != n; ++i)
        f[i] = M.mul(f[i], g[i]);
      intt(f, n);
      for (u32 i = 0; i != n; ++i)
        f[i] = M.get(f[i]);
    } else {
      const auto MX = _MX;
      for (u32 i = 0; i != n; i += 16) {
        MX.storeu(f + i, MX.trans(MX.loadu(f + i)));
        MX.storeu(g + i, MX.trans(MX.loadu(g + i)));
      }
      ntt(f, n), ntt(g, n);
      for (u32 i = 0; i != n; i += 16) {
        u32x16 fx = MX.loadu(f + i), gx = MX.loadu(g + i);
        MX.storeu(f + i, MX.mul(fx, gx));
      }
      intt(f, n);
      for (u32 i = 0; i != n; i += 16) {
        MX.storeu(f + i, MX.get(MX.loadu(f + i)));
      }
    }
  }
};

ALGO_END_NAMESPACE
