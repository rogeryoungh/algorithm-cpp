#ifndef ALGO_MATH_POLY_DIV10E_NTBLOCK
#define ALGO_MATH_POLY_DIV10E_NTBLOCK

#include "../../base.hpp"
#include "ntt.hpp"
#include "inv-10E-nt-block.hpp"
#include "nt-block-helper.hpp"

#include <algorithm>
#include <vector>
#include <iostream>

template <class ModT>
std::vector<ModT> poly_div_10E_block(std::span<const ModT> lhs, std::span<const ModT> rhs, u32 m) {
  if (lhs.empty() || rhs.empty())
    return {};
  if (m == 1)
    return {lhs[0] / rhs[0]};
  auto [n, u] = detail::nt_block_len(m);
  std::vector<ModT> x = poly_div_10E_block(lhs, rhs, n), h = poly_inv_10E_block(rhs, n);
  x.resize(n * u), h.resize(n * 2);
  std::vector<ModT> nf0(n * u * 2), ng0(n * u * 2);
  auto nf = detail::nt_block_split(nf0, n * 2);
  auto ng = detail::nt_block_split(ng0, n * 2);
  auto xk = detail::nt_block_split(x, n);

  ntt<ModT>(h);
  for (u32 k = 0; k < u; ++k) {
    std::copy(rhs.begin() + k * n, std::min(rhs.begin() + (k + 1) * n, rhs.end()), nf[k].begin());
    ntt<ModT>(nf[k]);
    if (k == 0)
      continue;
    std::copy(xk[k - 1].begin(), xk[k - 1].end(), ng[k - 1].begin());
    ntt<ModT>(ng[k - 1]);
    std::vector<ModT> psi(n * 2);
    for (u32 j = 0; j < k; ++j) {
      for (u32 i = 0; i < n; ++i)
        psi[i] -= (nf[k - j][i] + nf[k - 1 - j][i]) * ng[j][i];
      for (u32 i = n; i < n * 2; ++i)
        psi[i] -= (nf[k - j][i] - nf[k - 1 - j][i]) * ng[j][i];
    }
    intt<ModT>(psi);
    std::fill_n(psi.begin() + n, n, 0);
    for (u32 j = 0; j < std::min<u32>(n, lhs.size() - n * k); ++j)
      psi[j] += lhs[n * k + j];
    ntt<ModT>(psi);
    dot<ModT>(psi, h);
    intt<ModT>(psi);
    std::copy_n(psi.begin(), n, xk[k].begin());
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_DIV10E_NTBLOCK
