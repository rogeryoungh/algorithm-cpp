#ifndef ALGO_MATH_POLY_DOTS
#define ALGO_MATH_POLY_DOTS

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <span>

#ifndef ALGO_DISABLE_SIMD_AVX2

#include "../../other/modint/montgomery-x8.hpp"

template <class ModT>
static void dot_avx(std::span<ModT> f, std::span<const ModT> g) {
  using X8 = simd::M32x8<ModT>;
  using simd::i256::load, simd::i256::store;
  u32 n = f.size(), lf = u64(&f[0]) & 0x1f, i = 0;
  const auto loadg = lf == (u64(&g[0]) & 0x1f) ? load<true> : load<false>;
  for (; u64(&f[i]) & 0x1f; ++i) {
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

template <class ModT>
static void dot_v_avx(std::span<ModT> f, ModT v) {
  using X8 = simd::M32x8<ModT>;
  using simd::i256::load, simd::i256::store;
  u32 n = f.size(), i = 0;
  for (; u64(&f[i]) & 0x1f; ++i) {
    f[i] *= v;
  }
  X8 vx8 = X8::from(v);
  for (; i + 7 < n; i += 8) {
    X8 fi = load((simd::I256 *)&f[i]);
    fi *= vx8;
    store((simd::I256 *)&f[i], fi.v);
  }
  for (; i < n; ++i) {
    f[i] *= v;
  }
}

template <class ModT>
static void dot_avx(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  using X8 = simd::M32x8<ModT>;
  using simd::i256::load, simd::i256::store;
  u32 n = f.size(), lf = u64(&f[0]) & 0x1f, i = 0;
  const auto loadg = lf == (u64(&g[0]) & 0x1f) ? load<true> : load<false>;
  const auto storeg = lf == (u64(&dst[0]) & 0x1f) ? store<true> : store<false>;
  for (; u64(&f[i]) & 0x1f; ++i) {
    dst[i] = f[i] * g[i];
  }
  for (; i + 7 < n; i += 8) {
    X8 fi = load((simd::I256 *)&f[i]);
    X8 gi = loadg((simd::I256 *)&g[i]);
    fi *= gi;
    stored((simd::I256 *)&dst[i], fi.v);
  }
  for (; i < n; ++i) {
    dst[i] = f[i] * g[i];
  }
}

template <class ModT>
static void dot_v_avx(std::span<ModT> f, ModT v, std::span<ModT> dst) {
  using X8 = simd::M32x8<ModT>;
  using simd::i256::load, simd::i256::store;
  u32 n = f.size(), lf = u64(&f[0]) & 0x1f;
  const auto storeg = lf == (u64(&dst[0]) & 0x1f) ? store<true> : store<false>;
  u32 i = 0;
  for (; u64(&f[i]) & 0x1f; ++i) {
    dst[i] = f[i] * v;
  }
  X8 vx8 = X8::from(v);
  for (; i + 7 < n; i += 8) {
    X8 fi = load((simd::I256 *)&f[i]);
    fi *= vx8;
    stored((simd::I256 *)&dst[i], fi.v);
  }
  for (; i < n; ++i) {
    dst[i] = f[i] * v;
  }
}

#endif

template <class ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && f.size() > 16) {
    dot_avx(f, g, dst);
  } else {
#endif
    u32 n = dst.size();
    for (u32 i = 0; i < n; i++)
      dst[i] = f[i] * g[i];
#ifndef ALGO_DISABLE_SIMD_AVX2
  }
#endif
}

template <class ModT>
static void dot(std::span<ModT> f, std::span<const ModT> g) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && f.size() > 16) {
    dot_avx(f, g);
  } else {
#endif
    u32 n = f.size();
    for (u32 i = 0; i < n; i++)
      f[i] *= g[i];
#ifndef ALGO_DISABLE_SIMD_AVX2
  }
#endif
}

template <class ModT>
static void dot_v(std::span<ModT> f, ModT v, std::span<ModT> dst) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && f.size() > 16) {
    dot_v_avx(f, v, dst);
  } else {
#endif
    u32 n = dst.size();
    for (u32 i = 0; i < n; i++)
      dst[i] = f[i] * v;
#ifndef ALGO_DISABLE_SIMD_AVX2
  }
#endif
}

template <class ModT>
static void dot_v(std::span<ModT> f, ModT v) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && f.size() > 16) {
    dot_v_avx(f, v);
  } else {
#endif
    u32 n = f.size();
    for (u32 i = 0; i < n; i++)
      f[i] *= v;
#ifndef ALGO_DISABLE_SIMD_AVX2
  }
#endif
}

#endif // ALGO_MATH_POLY_DOTS
