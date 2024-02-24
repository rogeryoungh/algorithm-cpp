#pragma once

#include "../../number/montgomery.hpp"

#include <bit>
#include <array>

ALGO_BEGIN_NAMESPACE

template <class U>
struct NTTOriginalRadix2 {
  const Mont<U> _M;
  std::array<U, 64> rt{}, irt{}, rate2{}, irate2{};
  NTTOriginalRadix2(Mont<U> M, TI<U> G) : _M(std::move(M)) {
    const u32 rank2 = std::countr_zero(M.MOD - 1);
    rt[rank2] = M.qpow(M.trans(G), (M.MOD - 1) >> rank2);
    irt[rank2] = M.inv(rt[rank2]);
    for (u32 i = rank2; i != 0; --i) {
      rt[i - 1] = M.mul(rt[i], rt[i]);
      irt[i - 1] = M.mul(irt[i], irt[i]);
    }
    U prod = M.ONE, iprod = M.ONE;
    for (u32 i = 0; i != rank2 - 1; ++i) {
      rate2[i] = M.mul(prod, rt[i + 2]);
      irate2[i] = M.mul(iprod, irt[i + 2]);
      prod = M.mul(prod, irt[i + 2]);
      iprod = M.mul(iprod, rt[i + 2]);
    }
  }
  void ntt(U *f, usize n) {
    const auto M = _M;
    for (u32 l = n / 2; l != 0; l /= 2) {
      U r = M.ONE;
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          U x = f[i + j], y = M.mul(f[i + j + l], r);
          f[i + j] = M.add(x, y);
          f[i + j + l] = M.sub(x, y);
        }
        r = M.mul(r, rate2[std::countr_one(k)]);
      }
    }
  }
  void intt(U *f, usize n) {
    const auto M = _M;
    for (u32 l = 1; l != n; l *= 2) {
      U r = M.ONE;
      for (u32 i = 0, k = 0; i != n; i += l * 2, ++k) {
        for (u32 j = 0; j != l; ++j) {
          U x = f[i + j], y = f[i + j + l];
          f[i + j] = M.add(x, y);
          f[i + j + l] = M.mul(M.sub(x, y), r);
        }
        r = M.mul(r, irate2[std::countr_one(k)]);
      }
    }
    U ivn = M.trans(M.MOD - (M.MOD - 1) / n);
    for (u32 i = 0; i != n; ++i) {
      f[i] = M.mul(f[i], ivn);
    }
  }
  void conv(U *f, U *g, u32 n) {
    const auto M = _M;
    ntt(f, n);
    if (f != g)
      ntt(g, n);
    for (u32 i = 0; i != n; ++i)
      f[i] = M.mul(f[i], g[i]);
    intt(f, n);
  }
};

ALGO_END_NAMESPACE
