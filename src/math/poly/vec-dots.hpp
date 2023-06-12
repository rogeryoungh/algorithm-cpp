#ifndef ALGO_MATH_POLY_DOTS
#define ALGO_MATH_POLY_DOTS

#include "../../base.hpp"
#include "../../other/modint/modint-concept.hpp"

#include <span>

#ifndef ALGO_DISABLE_SIMD_AVX2
#include "../../other/modint/montgomery-x8.hpp"
#endif

template <class ModT, bool align = false, bool extra = true>
inline void vectorization_1(u32 n, ModT *f, auto &&op) {
  u32 i = 0;
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    using simd::i256::load, simd::i256::store;
    if constexpr (!align)
      for (; u64(&f[i]) & 0x1f; ++i)
        op(f[i]);
    for (; i + 7 < n; i += 8) {
      X8 fi = load((simd::I256 *)&f[i]);
      op(fi);
      store((simd::I256 *)&f[i], fi.v);
    }
  }
#endif
  if constexpr (extra)
    for (; i < n; ++i)
      op(f[i]);
}

template <class ModT, bool align = false, bool extra = true>
inline void vectorization_2(u32 n, ModT *f, const ModT *g, auto &&op) {
  u32 i = 0;
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    using simd::i256::load, simd::i256::store;
    if constexpr (!align)
      for (; u64(&f[i]) & 0x1f; ++i)
        op(f[i], g[i]);
    for (; i + 7 < n; i += 8) {
      X8 fi = load((simd::I256 *)&f[i]);
      X8 gi = load<align>((simd::I256 *)&g[i]);
      op(fi, gi);
      store((simd::I256 *)&f[i], fi.v);
    }
  }
#endif
  if constexpr (extra)
    for (; i < n; ++i)
      op(f[i], g[i]);
}

template <class ModT, bool align = false, bool extra = true>
void vectorization_3(u32 n, ModT *f, const ModT *g1, const ModT *g2, auto &&op) {
  u32 i = 0;
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    using simd::i256::load, simd::i256::store;
    if constexpr (!align)
      for (; u64(&f[i]) & 0x1f; ++i)
        op(f[i], g1[i], g2[i]);
    for (; i + 7 < n; i += 8) {
      X8 fi = load((simd::I256 *)&f[i]);
      X8 g1i = load<align>((simd::I256 *)&g1[i]);
      X8 g2i = load<align>((simd::I256 *)&g2[i]);
      op(fi, g1i, g2i);
      store((simd::I256 *)&f[i], fi.v);
    }
  }
#endif
  if constexpr (extra)
    for (; i < n; ++i)
      op(f[i], g1[i], g2[i]);
}

template <class ModT, bool align = false, bool extra = true>
void vectorization_4(u32 n, ModT *f, const ModT *g1, const ModT *g2, const ModT *g3, auto &&op) {
  u32 i = 0;
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    using simd::i256::load, simd::i256::store;
    if constexpr (!align)
      for (; u64(&f[i]) & 0x1f; ++i)
        op(f[i], g1[i], g2[i], g3[i]);
    for (; i + 7 < n; i += 8) {
      X8 fi = load((simd::I256 *)&f[i]);
      X8 g1i = load<align>((simd::I256 *)&g1[i]);
      X8 g2i = load<align>((simd::I256 *)&g2[i]);
      X8 g3i = load<align>((simd::I256 *)&g3[i]);
      op(fi, g1i, g2i, g3i);
      store((simd::I256 *)&f[i], fi.v);
    }
  }
#endif
  if constexpr (extra)
    for (; i < n; ++i)
      op(f[i], g1[i], g2[i], g3[i]);
}

template <class ModT, bool align = true, bool extra = true>
void dot(std::span<const ModT> f, std::span<const ModT> g, std::span<ModT> dst) {
  vectorization_3<ModT, align, extra>(f.size(), dst.data(), f.data(), g.data(), [](auto &di, auto fi, auto gi) {
    di = fi * gi;
  });
}

template <class ModT, bool align = true, bool extra = true>
void dot(std::span<ModT> f, std::span<const ModT> g) {
  vectorization_2<ModT, align, extra>(f.size(), f.data(), g.data(), [](auto &fi, auto gi) {
    fi *= gi;
  });
}

template <class ModT, bool align = true, bool extra = true>
static void dot_v(std::span<ModT> f, ModT v, std::span<ModT> dst) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  auto vx8 = simd::M32x8<ModT>::from(v);
  vectorization_2<ModT, align, extra>(dst.size(), dst.data(), f.data(), [v, vx8](auto &di, auto fi) {
    if constexpr (sizeof(fi) == sizeof(ModT))
      di = fi * v;
    else
      di = fi * vx8;
  });
#else
  vectorization_2<ModT, align, extra>(dst.size(), dst.data(), f.data(), [v](auto &di, auto fi) {
    di = fi * v;
  });
#endif
}

template <class ModT, bool align = true, bool extra = true>
static void dot_v(std::span<ModT> f, ModT v) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  auto vx8 = simd::M32x8<ModT>::from(v);
  vectorization_1<ModT, align, extra>(f.size(), f.data(), [v, vx8](auto &fi) {
    if constexpr (sizeof(fi) == sizeof(ModT))
      fi *= v;
    else
      fi *= vx8;
  });
#else
  vectorization_1<ModT, align, extra>(f.size(), f.data(), [v](auto &fi) {
    fi *= v;
  });
#endif
}

#endif // ALGO_MATH_POLY_DOTS
