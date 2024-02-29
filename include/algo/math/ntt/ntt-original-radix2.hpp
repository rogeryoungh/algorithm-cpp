#pragma once

#include "../../number/montgomery.hpp"
#include "./conv-ntt3.hpp"

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
    for (u32 l = n / 2; l >= 1; l /= 2) {
      U *fx = f, *fy = fx + l;
      for (u32 j = 0; j != l; ++j) {
        U x = fx[j], y = fy[j];
        fx[j] = M.add(x, y);
        fy[j] = M.sub(x, y);
      }
      U r = rate2[0];
      for (u32 i = l * 2, k = 1; i != n; i += l * 2, ++k) {
        fx = f + i, fy = fx + l;
        for (u32 j = 0; j != l; ++j) {
          U x = fx[j], y = M.mul(fy[j], r);
          fx[j] = M.add(x, y);
          fy[j] = M.sub(x, y);
        }
        r = M.mul(r, rate2[std::countr_one(k)]);
      }
    }
  }
  void intt(U *f, usize n) {
    const auto M = _M;
    U ivn = M.trans(M.MOD - (M.MOD - 1) / n);
    for (u32 l = 1; l <= n / 2; l *= 2) {
      U *fx = f, *fy = fx + l;
      for (u32 j = 0; j != l; ++j) {
        U x = fx[j], y = fy[j];
        if (l == n / 2)
          x = M.mul(x, ivn), y = M.mul(y, ivn); // div n here !!!
        fx[j] = M.add(x, y);
        fy[j] = M.sub(x, y);
      }
      U r = irate2[0];
      for (u32 i = l * 2, k = 1; i != n; i += l * 2, ++k) {
        fx = f + i, fy = fx + l;
        for (u32 j = 0; j != l; ++j) {
          U x = fx[j], y = fy[j];
          fx[j] = M.add(x, y);
          fy[j] = M.mul(M.sub(x, y), r);
        }
        r = M.mul(r, irate2[std::countr_one(k)]);
      }
    }
  }
  void conv(U *f, U *g, usize n) {
    conv_ntt3(f, g, n, *this);
  }
};

ALGO_END_NAMESPACE
