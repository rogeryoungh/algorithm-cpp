#pragma once

#include "../../number/montgomery.hpp"

#include <vector>

ALGO_BEGIN_NAMESPACE

template <class U>
struct NTTTwistedRadix2 {
  const Mont<U> _M;
  std::vector<U> rt, irt;
  NTTTwistedRadix2(Mont<U> M, TI<U> G) : _M(std::move(M)), rt(2), irt(2) {
    rt[0] = M.trans(G), irt[1] = rt[1] = M.ONE;
  }
  void prepare_root(usize m) {
    usize n = rt.size();
    if (n >= m)
      return;
    rt.resize(m), irt.resize(m);
    const auto M = _M;
    for (; n != m; n *= 2) {
      U w = M.qpow(rt[0], (M.MOD - 1) / n / 2), iw = M.inv(w);
      for (u32 i = n; i != n * 2; i += 2) {
        rt[i] = rt[i / 2];
        rt[i + 1] = M.mul(w, rt[i]);
        irt[i] = irt[i / 2];
        irt[i + 1] = M.mul(iw, irt[i]);
      }
    }
  }
  void ntt(U *f, usize n) {
    if (n <= 1)
      return;
    const auto M = _M;
    prepare_root(n);
    for (u32 l = n / 2; l >= 2; l /= 2) {
      U *w = rt.data() + l;
      for (u32 i = 0; i != n; i += l * 2) {
        U *fx = f + i, *fy = fx + l;
        for (u32 j = 0; j != l; ++j) {
          U x = fx[j], y = fy[j];
          fx[j] = M.add(x, y);
          fy[j] = M.mul(M.sub(x, y), w[j]);
        }
      }
    }
    for (u32 i = 0; i != n; i += 2) {
      U x = f[i + 0], y = f[i + 1];
      f[i + 0] = M.add(x, y);
      f[i + 1] = M.sub(x, y);
    }
  }
  void intt(U *f, usize n) {
    if (n <= 1)
      return;
    const auto M = _M;
    prepare_root(n);
    U ivn = M.trans(M.MOD - (M.MOD - 1) / n);
    for (u32 i = 0; i != n; i += 2) {
      U x = M.mul(f[i + 0], ivn), y = M.mul(f[i + 1], ivn);
      f[i + 0] = M.add(x, y);
      f[i + 1] = M.sub(x, y);
    }
    for (u32 l = 2; l <= n / 2; l *= 2) {
      U *w = irt.data() + l;
      for (u32 i = 0; i != n; i += l * 2) {
        U *fx = f + i, *fy = fx + l;
        for (u32 j = 0; j != l; ++j) {
          U x = fx[j], y = M.mul(fy[j], w[j]);
          fx[j] = M.add(x, y);
          fy[j] = M.sub(x, y);
        }
      }
    }
  }
  void conv(U *f, U *g, u32 n) {
    const auto M = _M;
    for (u32 i = 0; i != n; ++i)
      f[i] = M.trans(f[i]), g[i] = M.trans(g[i]);
    ntt(f, n), ntt(g, n);
    for (u32 i = 0; i != n; ++i)
      f[i] = M.mul(f[i], g[i]);
    intt(f, n);
    for (u32 i = 0; i != n; ++i)
      f[i] = M.get(f[i]);
  }
};

ALGO_END_NAMESPACE
