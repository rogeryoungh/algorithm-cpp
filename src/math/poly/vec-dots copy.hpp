#ifndef ALGO_MATH_POLY_DOTS
#define ALGO_MATH_POLY_DOTS

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <span>

template <static_modint_concept ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  u32 n = dst.size();
  for (u32 i = 0; i < n; i++)
    dst[i] = f[i] * g[i];
}

template <static_modint_concept ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g) {
  u32 n = f.size();
  for (u32 i = 0; i < n; i++)
    f[i] *= g[i];
}

#endif // ALGO_MATH_POLY_DOTS
