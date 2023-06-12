#ifndef ALGO_MATH_POLY_INV10E_NTBLOCK
#define ALGO_MATH_POLY_INV10E_NTBLOCK

#include "poly-def.hpp"
#include "nt-block-helper.hpp"

template <class ModT>
AVec<ModT> poly_inv_10E_block(std::span<const ModT> self, u32 m) {
  if (m == 1)
    return {self[0].inv()};
  auto [n, u] = detail::nt_block_len(m);
  AVec<ModT> x = poly_inv_10E_block(self, n);
  x.resize(n * u);
  AVec<ModT> nf0(n * u * 2), ng0(n * u * 2);
  auto nf = detail::nt_block_split(nf0, n * 2);
  auto ng = detail::nt_block_split(ng0, n * 2);
  auto xk = detail::nt_block_split(x, n);

  for (u32 k = 0; k < u; ++k) {
    std::copy(self.begin() + k * n, std::min(self.begin() + (k + 1) * n, self.end()), nf[k].begin());
    ntt<ModT>(nf[k]);
    if (k == 0)
      continue;
    std::copy(xk[k - 1].begin(), xk[k - 1].end(), ng[k - 1].begin());
    ntt<ModT>(ng[k - 1]);
    AVec<ModT> psi(n * 2);
    for (u32 j = 0; j < k; ++j) {
      auto psi_p = psi.data(), nf1_p = nf[k - j].data(), nf2_p = nf[k - 1 - j].data(), ng_p = ng[j].data();
      const auto fn1 = []<class T>(T &pi, T nf1i, T nf2i, T ngi) {
        pi -= T::addmul(nf1i, nf2i, ngi);
      };
      vectorization_4<ModT, true>(n, psi_p, nf1_p, nf2_p, ng_p, fn1);
      const auto fn2 = []<class T>(T &pi, T nf1i, T nf2i, T ngi) {
        pi -= T::submul(nf1i, nf2i, ngi);
      };
      vectorization_4<ModT, true>(n, psi_p + n, nf1_p + n, nf2_p + n, ng_p + n, fn2);
    }
    intt<ModT>(psi);
    std::fill_n(psi.begin() + n, n, 0);
    ntt<ModT>(psi);
    dot<ModT>(psi, ng[0]);
    intt<ModT>(psi);
    std::copy_n(psi.begin(), n, xk[k].begin());
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INV10E_NTBLOCK
