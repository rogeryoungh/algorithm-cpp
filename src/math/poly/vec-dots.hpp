#ifndef ALGO_MATH_POLY_DOTS
#define ALGO_MATH_POLY_DOTS

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <span>

template <class ModT>
static void dot_basic(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  u32 n = dst.size();
  for (u32 i = 0; i < n; i++)
    dst[i] = f[i] * g[i];
}

template <class ModT>
static void dot_basic(std::span<ModT> f, std::span<const ModT> g) {
  u32 n = f.size();
  for (u32 i = 0; i < n; i++)
    f[i] *= g[i];
}

#ifndef ALGO_DISABLE_SIMD_AVX2

#include "../../other/modint/montgomery-x8.hpp"

template <class ModT>
static void dot_avx(std::span<ModT> f, std::span<const ModT> g) {
  u32 n8 = f.size();
  u32 i = 0;
  using X8 = simd::M32x8<ModT>;
  for (; i + 7 < n8; i += 8) {
    X8 fi = X8::load((simd::I256 *)&f[i]);
    X8 gi = X8::load((simd::I256 *)&g[i]);
    fi *= gi;
    fi.store((simd::I256 *)&f[i]);
  }
  for (; i < n8; i++)
    f[i] *= g[i];
}

template <class ModT>
static void dot_avx(std::span<simd::I256> f, std::span<const ModT> g, std::span<ModT> dst) {
  u32 n = dst.size();
  u32 i = 0;
  using X8 = simd::M32x8<ModT>;
  for (; i + 7 < n; i += 8) {
    X8 fi = X8::load((simd::I256 *)&f[i]);
    X8 gi = X8::load((simd::I256 *)&g[i]);
    X8 di = fi * gi;
    di.store((simd::I256 *)&dst[i]);
  }
  for (; i < n; i++)
    dst[i] = f[i] * g[i];
}

template <class ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  if constexpr (montgomery_modint_concept<ModT>) {
    dot_avx(f, g, dst);
  } else {
    dot_basic(f, g, dst);
  }
}

template <class ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g) {
  if constexpr (montgomery_modint_concept<ModT>) {
    dot_avx(f, g);
  } else {
    dot_basic(f, g);
  }
}

#else

template <class ModT>
void dot(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  dot_basic(f, g, dst);
}

template <class ModT>
void dot(std::span<ModT> f, std::span<const ModT> g) {
  dot_basic(f, g);
}

#endif

#endif // ALGO_MATH_POLY_DOTS
