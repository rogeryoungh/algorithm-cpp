#ifndef ALGO_MATH_POLY_EXP14E_NTBLOCK
#define ALGO_MATH_POLY_EXP14E_NTBLOCK

#include "../../base.hpp"
#include "ntt.hpp"
#include "vec-dots.hpp"
#include "nt-block-helper.hpp"
#include "inv-10E-nt-block.hpp"
#include "../constant/prepare-inv.hpp"

#include <algorithm>
#include <vector>
#include <iostream>

template <static_modint_concept ModT>
std::vector<ModT> poly_exp_14E_block(std::span<const ModT> self, u32 m) {
  if (m == 1)
    return {1};
  auto [n, u] = detail::nt_block_len(m);
  std::vector<ModT> x = poly_exp_14E_block(self, n), h = poly_inv_10E_block<ModT>(x, n);
  x.resize(n * u), h.resize(n * 2);
  std::vector<ModT> nf0(n * u * 2), ng0(n * u * 2);
  auto &iv = prepare_inv<ModT>(n * u);
  auto nf = detail::nt_block_split(nf0, n * 2);
  auto ng = detail::nt_block_split(ng0, n * 2);
  auto xk = detail::nt_block_split(x, n);

  ntt<ModT>(h);
  for (u32 k = 0; k < u; ++k) {
    std::copy(self.begin() + k * n, std::min(self.begin() + (k + 1) * n, self.end()), nf[k].begin());
    for (u32 i = 0; i < n; ++i)
      nf[k][i] *= k * n + i;
    ntt<ModT>(nf[k]);
    if (k == 0)
      continue;
    std::copy(xk[k - 1].begin(), xk[k - 1].end(), ng[k - 1].begin());
    ntt<ModT>(ng[k - 1]);
    std::vector<ModT> psi(n * 2);
    for (u32 j = 0; j < k; ++j) {
      for (u32 i = 0; i < n; ++i)
        psi[i] += (nf[k - j][i] + nf[k - 1 - j][i]) * ng[j][i];
      for (u32 i = n; i < n * 2; ++i)
        psi[i] += (nf[k - j][i] - nf[k - 1 - j][i]) * ng[j][i];
    }
    intt<ModT>(psi);
    std::fill_n(psi.begin() + n, n, 0);
    ntt<ModT>(psi);
    dot<ModT>(psi, h);
    intt<ModT>(psi);
    std::fill_n(psi.begin() + n, n, 0);
    dot<ModT>({psi.begin(), n}, {iv.begin() + n * k, n});
    ntt<ModT>(psi);
    dot<ModT>(psi, ng[0]);
    intt<ModT>(psi);
    std::copy_n(psi.begin(), n, xk[k].begin());
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_EXP14E_NTBLOCK
