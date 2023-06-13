#ifndef ALGO_MATH_POLY_DIV10E_NTBLOCK
#define ALGO_MATH_POLY_DIV10E_NTBLOCK

#include "nt-block-helper.hpp"
#include "poly-def.hpp"

template <class ModT, auto poly_inv>
AVec<ModT> poly_div_10E_block(std::span<const ModT> lhs, std::span<const ModT> rhs, u32 m) {
  if (lhs.empty() || rhs.empty())
    return {};
  if (m == 1)
    return {lhs[0] / rhs[0]};
  auto [n, u] = detail::nt_block_len(m);
  AVec<ModT> x = poly_div_10E_block<ModT, poly_inv>(lhs, rhs, n), h = poly_inv(rhs, n);
  x.resize(n * u), h.resize(n * 2);
  AVec<ModT> nf0(n * u * 2), ng0(n * u * 2);
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
    AVec<ModT> psi(n * 2);
    for (u32 j = 0; j < k; ++j) {
      auto psi_p = psi.data(), nf1_p = nf[k - j].data(), nf2_p = nf[k - 1 - j].data(), ng_p = ng[j].data();
      const auto fn1 = []<class T>(T &pi, T nf1i, T nf2i, T ngi) {
        pi -= (nf1i + nf2i) * ngi;
      };
      vectorization_4<ModT, true>(n, psi_p, nf1_p, nf2_p, ng_p, fn1);
      const auto fn2 = []<class T>(T &pi, T nf1i, T nf2i, T ngi) {
        pi -= (nf1i - nf2i) * ngi;
      };
      vectorization_4<ModT, true>(n, psi_p + n, nf1_p + n, nf2_p + n, ng_p + n, fn2);
    }
    intt<ModT>(psi);
    std::fill_n(psi.begin() + n, n, 0);
    vectorization_2<ModT, true>(
        std::min<u32>(n, lhs.size() - n * k), psi.data(), lhs.data() + n * k, []<class T>(T &pi, T li) {
          pi += li;
        });
    ntt<ModT>(psi);
    dot<ModT, true>(psi, h);
    intt<ModT>(psi);
    std::copy_n(psi.begin(), n, xk[k].begin());
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_DIV10E_NTBLOCK
