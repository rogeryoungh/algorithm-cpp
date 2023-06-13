#ifndef ALGO_MATH_POLY_SQRT8E_NTBLOCK
#define ALGO_MATH_POLY_SQRT8E_NTBLOCK

#include "poly-def.hpp"
#include "nt-block-helper.hpp"

template <class ModT, auto poly_inv>
AVec<ModT> poly_sqrt_8E_block(std::span<const ModT> self, u32 m, const ModT &x0) {
  if (m == 1)
    return {x0};
  auto [n, u] = detail::nt_block_len(m);
  AVec<ModT> x = poly_sqrt_8E_block<ModT, poly_inv>(self, n, x0), h = poly_inv(x, n);
  x.resize(n * u), h.resize(n * 2);
  AVec<ModT> ng0(n * u * 2);
  auto ng = detail::nt_block_split(ng0, n * 2);
  auto xk = detail::nt_block_split(x, n);

  ntt<ModT>(h);
  for (u32 k = 1; k < u; ++k) {
    std::copy(xk[k - 1].begin(), xk[k - 1].end(), ng[k - 1].begin());
    ntt<ModT>(ng[k - 1]);
    AVec<ModT> psi(n * 2);
    for (u32 j = 0; j < k; ++j) {
      auto psi_p = psi.data(), ng1_p = ng[k - j].data(), ng2_p = ng[k - 1 - j].data(), ng_p = ng[j].data();
      if (j == 0) {
        const auto fn1 = []<class T>(T &pi, T ng2i, T ngi) {
          pi += ng2i * ngi;
        };
        vectorization_3<ModT, true>(n, psi_p, ng2_p, ng_p, fn1);
        const auto fn2 = []<class T>(T &pi, T ng2i, T ngi) {
          pi -= ng2i * ngi;
        };
        vectorization_3<ModT, true>(n, psi_p + n, ng2_p + n, ng_p + n, fn2);

      } else {
        const auto fn1 = []<class T>(T &pi, T ng1i, T ng2i, T ngi) {
          pi += (ng1i + ng2i) * ngi;
        };
        vectorization_4<ModT, true>(n, psi_p, ng1_p, ng2_p, ng_p, fn1);
        const auto fn2 = []<class T>(T &pi, T ng1i, T ng2i, T ngi) {
          pi += (ng1i - ng2i) * ngi;
        };
        vectorization_4<ModT, true>(n, psi_p + n, ng1_p + n, ng2_p + n, ng_p + n, fn2);
      }
    }
    intt<ModT>(psi);
    std::fill_n(psi.begin() + n, n, 0);
    vectorization_2(std::min<u32>(n, self.size() - n * k), psi.data(), self.data() + n * k, []<class T>(T &pi, T si) {
      pi = (si - pi).shift2();
    });
    ntt<ModT>(psi);
    dot<ModT>(psi, h);
    intt<ModT>(psi);
    std::copy_n(psi.begin(), n, xk[k].begin());
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_SQRT8E_NTBLOCK
