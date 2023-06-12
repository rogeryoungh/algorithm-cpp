#ifndef ALGO_MATH_POLY_SQRT8E_NTBLOCK
#define ALGO_MATH_POLY_SQRT8E_NTBLOCK

#include "poly-def.hpp"
#include "nt-block-helper.hpp"

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

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
#ifndef ALGO_DISABLE_SIMD_AVX2
    if (montgomery_modint_concept<ModT> && n > 16) {
      using X8 = simd::M32x8<ModT>;
      auto *psix8 = reinterpret_cast<X8 *>(psi.data());
      u32 nx8 = n / 8;
      for (u32 j = 0; j < k; ++j) {
        auto *ng1x8 = reinterpret_cast<X8 *>(ng[k - j].data());
        auto *ng2x8 = reinterpret_cast<X8 *>(ng[k - j - 1].data());
        auto *ngx8 = reinterpret_cast<X8 *>(ng[j].data());
        if (j == 0) {
          for (u32 i = 0; i < nx8; ++i)
            psix8[i] += ng2x8[i] * ngx8[i];
          for (u32 i = nx8; i < nx8 * 2; ++i)
            psix8[i] -= ng2x8[i] * ngx8[i];
        } else {
          for (u32 i = 0; i < nx8; ++i)
            psix8[i] += (ng1x8[i] + ng2x8[i]) * ngx8[i];
          for (u32 i = nx8; i < nx8 * 2; ++i)
            psix8[i] += (ng1x8[i] - ng2x8[i]) * ngx8[i];
        }
      }
    } else {
#endif
      for (u32 j = 0; j < k; ++j) {
        if (j == 0) {
          for (u32 i = 0; i < n; i++)
            psi[i] += ng[k - 1 - j][i] * ng[j][i];
          for (u32 i = n; i < n * 2; i++)
            psi[i] -= ng[k - 1 - j][i] * ng[j][i];
        } else {
          for (u32 i = 0; i < n; i++)
            psi[i] += (ng[k - j][i] + ng[k - 1 - j][i]) * ng[j][i];
          for (u32 i = n; i < n * 2; i++)
            psi[i] += (ng[k - j][i] - ng[k - 1 - j][i]) * ng[j][i];
        }
      }
#ifndef ALGO_DISABLE_SIMD_AVX2
    }
#endif
    intt<ModT>(psi);
    std::fill_n(psi.begin() + n, n, 0);
    for (u32 j = 0; j < std::min<u32>(n, self.size() - n * k); j++)
      psi[j] = (self[n * k + j] - psi[j]).shift2();
    ntt<ModT>(psi);
    dot<ModT>(psi, h);
    intt<ModT>(psi);
    std::copy_n(psi.begin(), n, xk[k].begin());
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_SQRT8E_NTBLOCK
