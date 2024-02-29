#pragma once

#include "../../base.hpp"

ALGO_BEGIN_NAMESPACE

template <class U, class NTT>
void conv_ntt3(U *f, U *g, usize n, NTT &ntt) {
  const auto M = ntt._M;
  for (u32 i = 0; i != n; ++i)
    f[i] = M.trans(f[i]);
  if (f != g)
    for (u32 i = 0; i != n; ++i)
      g[i] = M.trans(g[i]);
  ntt.ntt(f, n);
  if (f != g)
    ntt.ntt(g, n);
  for (u32 i = 0; i != n; ++i)
    f[i] = M.mul(f[i], g[i]);
  ntt.intt(f, n);
  for (u32 i = 0; i != n; ++i)
    f[i] = M.get(f[i]);
}

ALGO_END_NAMESPACE
