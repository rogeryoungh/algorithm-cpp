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
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    if constexpr (!align)
      for (; u64(f) & 0x1f; --n)
        op(*f++);
    auto fx8 = std::bit_cast<X8 *>(f);
    for (; n >= 8; n -= 8) {
      op(*fx8++);
    }
    f = std::bit_cast<ModT *>(fx8);
  }
#endif
  if constexpr (extra)
    for (; n > 0; --n)
      op(*f++);
}

template <class ModT, bool align = false, bool extra = true>
inline void vectorization_2(u32 n, ModT *f, const ModT *g1, auto &&op) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    using i256a = simd::i256a<align>;
    if constexpr (!align)
      for (; u64(f) & 0x1f; --n)
        op(*f++, *g1++);
    auto fx8 = std::bit_cast<X8 *>(f);
    auto g1x8 = std::bit_cast<const i256a *>(g1);
    for (; n >= 8; n -= 8) {
      op(*fx8++, X8(*g1x8++));
    }
    f = std::bit_cast<ModT *>(fx8);
    g1 = std::bit_cast<const ModT *>(g1x8);
  }
#endif
  if constexpr (extra)
    for (; n > 0; --n)
      op(*f++, *g1++);
}

template <class ModT, bool align = false, bool extra = true>
void vectorization_3(u32 n, ModT *f, const ModT *g1, const ModT *g2, auto &&op) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    using i256a = simd::i256a<align>;
    if constexpr (!align)
      for (; u64(&f) & 0x1f; n -= 8)
        op(*f++, *g1++, *g2++);
    auto fx8 = std::bit_cast<X8 *>(f);
    auto g1x8 = std::bit_cast<const i256a *>(g1);
    auto g2x8 = std::bit_cast<const i256a *>(g2);
    for (; n >= 8; n -= 8) {
      op(*fx8++, X8(*g1x8++), X8(*g2x8++));
    }
    f = std::bit_cast<ModT *>(fx8);
    g1 = std::bit_cast<const ModT *>(g1x8);
    g2 = std::bit_cast<const ModT *>(g2x8);
  }
#endif
  if constexpr (extra)
    for (; n > 0; --n)
      op(*f++, *g1++, *g2++);
}

template <class ModT, bool align = false, bool extra = true>
void vectorization_4(u32 n, ModT *f, const ModT *g1, const ModT *g2, const ModT *g3, auto &&op) {
#ifndef ALGO_DISABLE_SIMD_AVX2
  if (montgomery_modint_concept<ModT> && n > 16) {
    using X8 = simd::M32x8<ModT>;
    using i256a = simd::i256a<align>;
    if constexpr (!align)
      for (; u64(&f) & 0x1f; n -= 8)
        op(*f++, *g1++, *g2++, *g3++);
    auto fx8 = std::bit_cast<X8 *>(f);
    auto g1x8 = std::bit_cast<const i256a *>(g1);
    auto g2x8 = std::bit_cast<const i256a *>(g2);
    auto g3x8 = std::bit_cast<const i256a *>(g3);
    for (; n >= 8; n -= 8) {
      op(*fx8++, X8(*g1x8++), X8(*g2x8++), X8(*g3x8++));
    }
    f = std::bit_cast<ModT *>(fx8);
    g1 = std::bit_cast<const ModT *>(g1x8);
    g2 = std::bit_cast<const ModT *>(g2x8);
    g3 = std::bit_cast<const ModT *>(g3x8);
  }
#endif
  if constexpr (extra)
    for (; n > 0; --n)
      op(*f++, *g1++, *g2++, *g3++);
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
