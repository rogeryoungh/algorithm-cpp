#ifndef ALGO_MATH_POLY_DOTS
#define ALGO_MATH_POLY_DOTS

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"
#include "../../other/avx.hpp"

#include <span>

// template <static_modint_concept ModT>
// static void dot(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
//   u32 n = dst.size();
//   for (u32 i = 0; i < n; i++)
//     dst[i] = f[i] * g[i];
// }

// template <static_modint_concept ModT>
// static void dot(std::span<ModT> f, std::span<const ModT> g) {
//   u32 n = f.size();
//   for (u32 i = 0; i < n; i++)
//     f[i] *= g[i];
// }

template <montgomery_modint_concept ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g) {
  u32 n = f.size();
  u32 i = 0;
  for (; i + 7 < n; i += 8) {
    auto pf = (u8x32 *)&f[i];
    auto pg = (const u8x32 *)&g[i];
    auto di = m32x8(pf) * m32x8(pg);
    di.v.store(pf);
  }
  for (; i < n; i++)
    f[i] *= g[i];
}

template <montgomery_modint_concept ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  u32 n = dst.size();
  u32 i = 0;
  for (; i + 7 < n; i += 8) {
    auto pf = (u8x32 *)&f[i];
    auto pg = (const u8x32 *)&g[i];
    auto pd = (u8x32 *)&dst[i];
    auto di = m32x8(pf) * m32x8(pg);
    di.v.store(pd);
  }
  for (; i < n; i++)
    dst[i] = f[i] * g[i];
}

#endif // ALGO_MATH_POLY_DOTS
