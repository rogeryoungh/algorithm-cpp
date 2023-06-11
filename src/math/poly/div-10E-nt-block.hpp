#ifndef ALGO_MATH_POLY_DIV10E_NTBLOCK
#define ALGO_MATH_POLY_DIV10E_NTBLOCK

#include "poly-def.hpp"
#include "inv-10E-nt-block.hpp"
#include "nt-block-helper.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

template <class ModT>
AVec<ModT> poly_div_10E_block(std::span<const ModT> lhs, std::span<const ModT> rhs, u32 m) {
  if (lhs.empty() || rhs.empty())
    return {};
  if (m == 1)
    return {lhs[0] / rhs[0]};
  auto [n, u] = detail::nt_block_len(m);
  AVec<ModT> x = poly_div_10E_block(lhs, rhs, n), h = poly_inv_10E_block(rhs, n);
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
