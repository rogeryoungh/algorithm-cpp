#ifndef ALGO_H_MATH_AVX2_NTT_BASE
#define ALGO_H_MATH_AVX2_NTT_BASE

#include "../../base.hpp"
#include "../../number/avx2/mont32x8.hpp"
#include "../../number/mont-vec-dots.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTT32BaseAVX2 {
  void dot(ModT *f, const ModT *g, u32 n) {
    if (n < 8) {
      mont_dot(f, g, n);
    } else {
      M32x8<ModT>::dot(f, g, n);
    }
  }
  void rescale(ModT *f, u32 n) {
    ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
    if (n < 8) {
      mont_dot1(f, n, ivn);
    } else {
      M32x8<ModT>::dot1(f, n, ivn);
    }
  }
};

ALGO_END_NAMESPACE

#endif
