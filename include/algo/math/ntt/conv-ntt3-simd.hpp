#pragma once

#include "../../base.hpp"

ALGO_BEGIN_NAMESPACE

template <u32 vsize, class U, class NTT>
void conv_ntt3_simd(U *f, U *g, usize n, NTT &ntt) {
  const auto MX = ntt._MX;
  for (u32 i = 0; i != n; i += vsize)
    MX.store(f + i, MX.trans(MX.load(f + i)));
  if (f != g)
    for (u32 i = 0; i != n; i += vsize)
      MX.store(g + i, MX.trans(MX.load(g + i)));
  ntt.ntt(f, n);
  if (f != g)
    ntt.ntt(g, n);
  for (u32 i = 0; i != n; i += vsize)
    MX.store(f + i, MX.mul(MX.load(f + i), MX.load(g + i)));
  ntt.intt(f, n);
  for (u32 i = 0; i != n; i += vsize)
    MX.store(f + i, MX.get(MX.load(f + i)));
}

ALGO_END_NAMESPACE
