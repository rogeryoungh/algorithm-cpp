#pragma once

#include "../../number/mont32x8.hpp"
#include <array>
#include <bit>

ALGO_BEGIN_NAMESPACE

struct NTT32OriginalRadix2AVX2 {
  std::array<u32, 32> rt, irt, rate2, irate2;
  u32x8 rate4ix[32], irate4ix[32];
  u32x8 _rt2, _irt2, _rt4, _irt4;
  Mont32x8 _MX;
  NTT32OriginalRadix2AVX2(Mont32 M, u32 G) : _MX(M) {
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
    u32 arr[8];
    auto rotate = [&M, &arr](u32 x) {
      for (u32 i = 0; i != 8; ++i)
        arr[i] = i == 0 ? M.ONE : M.mul(x, arr[i - 1]);
    };
    for (u32 i = 0; i != rank2 - 3; ++i) {
      rotate(M.mul(prod, rt[i + 4]));
      rate4ix[i] = _MX.loadu(arr);
      rotate(M.mul(iprod, irt[i + 4]));
      irate4ix[i] = _MX.loadu(arr);
      prod = M.mul(prod, irt[i + 4]);
      iprod = M.mul(iprod, rt[i + 4]);
    }
    auto rotatex = [&M, &arr](u32 x, u32 k) {
      for (u32 i = 0; i != 8; i += k)
        for (u32 j = 0; j != k; ++j)
          arr[i + j] = (j <= k / 2) ? M.ONE : M.mul(x, arr[i + j - 1]);
    };
    rotatex(rt[2], 4), _rt2 = _MX.loadu(arr);
    rotatex(irt[2], 4), _irt2 = _MX.loadu(arr);
    rotatex(rt[3], 8), _rt4 = _MX.loadu(arr);
    rotatex(irt[3], 8), _irt4 = _MX.loadu(arr);
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
    if (n < 8)
      return ntt_small(f, n);
    const auto MX = _MX;
    const auto M = MX.M;
    for (u32 l = n / 2; l >= 8; l /= 2) {
      u32 *f0 = f, *f1 = f + l;
      for (u32 j = 0; j != l; j += 8) {
        u32x8 x = MX.loadu(f0 + j), y = MX.loadu(f1 + j);
        MX.storeu(f0 + j, MX.add(x, y));
        MX.storeu(f1 + j, MX.sub(x, y));
      }
      u32 r = rate2[0];
      for (u32 i = l * 2, k = 1; i != n; i += l * 2, ++k) {
        u32x8 rx = MX.set1(r);
        f0 = f + i, f1 = f0 + l;
        for (u32 j = 0; j != l; j += 8) {
          u32x8 x = MX.loadu(f0 + j), y = MX.mul(rx, MX.loadu(f1 + j));
          MX.storeu(f0 + j, MX.add(x, y));
          MX.storeu(f1 + j, MX.sub(x, y));
        }
        r = M.mul(r, rate2[std::countr_one(k)]);
      }
    }
    u32x8 rtix = MX.ONE, rt2 = _rt2, rt4 = _rt4;
    for (u32 i = 0; i != n; i += 8) {
      u32x8 fi = MX.mul(rtix, MX.loadu(f + i)), a, b;
      a = MX.neg<0xf0>(fi), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      fi = MX.mul(rt4, MX.add(a, b));
      a = MX.neg<0xcc>(fi), b = u32x8_shuffle<0x4e>(fi);
      fi = MX.mul(rt2, MX.add(a, b));
      a = MX.neg<0xaa>(fi), b = u32x8_shuffle<0xb1>(fi);
      MX.storeu(f + i, MX.add(a, b));
      rtix = MX.mul(rtix, rate4ix[std::countr_one(i / 8)]);
    }
  }
  void intt(u32 *f, usize n) {
    if (n < 8)
      return intt_small(f, n);
    const auto MX = _MX;
    const auto M = MX.M;
    u32x8 rtix = MX.set1(M.trans(M.MOD - (M.MOD - 1) / n));
    u32x8 irt2 = _irt2, irt4 = _irt4;
    for (u32 i = 0; i != n; i += 8) {
      u32x8 fi = MX.loadu(f + i), a, b;
      a = MX.neg<0xaa>(fi), b = u32x8_shuffle<0xb1>(fi);
      fi = MX.mul(irt2, MX.add(a, b));
      a = MX.neg<0xcc>(fi), b = u32x8_shuffle<0x4e>(fi);
      fi = MX.mul(irt4, MX.add(a, b));
      a = MX.neg<0xf0>(fi), b = _mm256_permute2x128_si256(fi, fi, 0b01);
      MX.storeu(f + i, MX.mul(MX.add(a, b), rtix));
      rtix = MX.mul(rtix, irate4ix[std::countr_one(i / 8)]);
    }
    for (u32 l = 8; l <= n / 2; l *= 2) {
      u32 *f0 = f, *f1 = f + l;
      for (u32 j = 0; j != l; j += 8) {
        u32x8 x = MX.loadu(f0 + j), y = MX.loadu(f1 + j);
        MX.storeu(f0 + j, MX.add(x, y));
        MX.storeu(f1 + j, MX.sub(x, y));
      }
      u32 r = irate2[0];
      for (u32 i = l * 2, k = 1; i != n; i += l * 2, ++k) {
        u32x8 rx = MX.set1(r);
        f0 = f + i, f1 = f0 + l;
        for (u32 j = 0; j != l; j += 8) {
          u32x8 x = MX.loadu(f0 + j), y = MX.loadu(f1 + j);
          MX.storeu(f0 + j, MX.add(x, y));
          MX.storeu(f1 + j, MX.mul(MX.sub(x, y), rx));
        }
        r = M.mul(r, irate2[std::countr_one(k)]);
      }
    }
  }
  void conv(u32 *f, u32 *g, u32 n) {
    if (n < 8) {
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
      for (u32 i = 0; i != n; i += 8) {
        MX.storeu(f + i, MX.trans(MX.loadu(f + i)));
        MX.storeu(g + i, MX.trans(MX.loadu(g + i)));
      }
      ntt(f, n), ntt(g, n);
      for (u32 i = 0; i != n; i += 8) {
        u32x8 fx = MX.loadu(f + i), gx = MX.loadu(g + i);
        MX.storeu(f + i, MX.mul(fx, gx));
      }
      intt(f, n);
      for (u32 i = 0; i != n; i += 8) {
        MX.storeu(f + i, MX.get(MX.loadu(f + i)));
      }
    }
  }
};

ALGO_END_NAMESPACE
