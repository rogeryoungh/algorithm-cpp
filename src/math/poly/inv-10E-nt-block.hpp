#ifndef ALGO_MATH_POLY_INV10E_NTBLOCK
#define ALGO_MATH_POLY_INV10E_NTBLOCK

#include "poly-def.hpp"
#include "nt-block-helper.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

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
#ifndef ALGO_DISABLE_SIMD_AVX2
    if (montgomery_modint_concept<ModT> && n > 16) {
      using X8 = simd::M32x8<ModT>;
      auto *psix8 = reinterpret_cast<X8 *>(psi.data());
      u32 nx8 = n / 8;
      for (u32 j = 0; j < k; ++j) {
        auto *nf1x8 = reinterpret_cast<X8 *>(nf[k - j].data());
        auto *nf2x8 = reinterpret_cast<X8 *>(nf[k - j - 1].data());
        auto *ngx8 = reinterpret_cast<X8 *>(ng[j].data());
        for (u32 i = 0; i < nx8; ++i)
          psix8[i] -= (nf1x8[i] + nf2x8[i]) * ngx8[i];
        for (u32 i = nx8; i < nx8 * 2; ++i)
          psix8[i] -= (nf1x8[i] - nf2x8[i]) * ngx8[i];
      }
    } else {
#endif
      for (u32 j = 0; j < k; ++j) {
        for (u32 i = 0; i < n; ++i)
          psi[i] -= (nf[k - j][i] + nf[k - 1 - j][i]) * ng[j][i];
        for (u32 i = n; i < n * 2; ++i)
          psi[i] -= (nf[k - j][i] - nf[k - 1 - j][i]) * ng[j][i];
      }
#ifndef ALGO_DISABLE_SIMD_AVX2
    }
#endif
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
