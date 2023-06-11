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
  using X8 = simd::M32x8<ModT>;
  using simd::i256::load, simd::i256::store;
  u32 n = f.size(), lf = u64(f.data()) & 0x1f;
  const auto loadg = lf == (u64(g.data()) & 0x1f) ? load<true> : load<false>;
  if (n < 16) {
    dot_basic(f, g);
  } else {
    u32 i = 0;
    for (; i < lf; ++i) {
      f[i] *= g[i];
    }
    for (; i + 7 < n; i += 8) {
      X8 fi = load((simd::I256 *)&f[i]);
      X8 gi = loadg((simd::I256 *)&g[i]);
      fi *= gi;
      store((simd::I256 *)&f[i], fi.v);
    }
    for (; i < n; ++i) {
      f[i] *= g[i];
    }
  }
}

template <class ModT>
static void dot_avx(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  using X8 = simd::M32x8<ModT>;
  using simd::i256::load, simd::i256::store;
  u32 n = f.size(), lf = u64(f.data()) & 0x1f;
  const auto loadg = lf == (u64(g.data()) & 0x1f) ? load<true> : load<false>;
  const auto storeg = lf == (u64(dst.data()) & 0x1f) ? store<true> : store<false>;
  if (n < 16) {
    dot_basic(f, g);
  } else {
    u32 i = 0;
    for (; i < lf; ++i) {
      dst[i] = f[i] * g[i];
    }
    for (; i + 7 < n; i += 8) {
      X8 fi = load((simd::I256 *)&f[i]);
      X8 gi = loadg((simd::I256 *)&g[i]);
      fi *= gi;
      stored((simd::I256 *)&dst[i], fi.v);
    }
    for (; i < n; ++i) {
      f[i] *= g[i];
    }
  }
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
