#ifndef ALGO_H_MATH_NTT_BASE
#define ALGO_H_MATH_NTT_BASE

#include "../base.hpp"
#include "../modular/mont-vec-dots.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
struct NTTBase {
  void dot(ModT *f, const ModT *g, u32 n) {
    mont_dot(f, g, n);
  }
  void rescale(ModT *f, u32 n) {
    ModT ivn = ModT::MOD - (ModT::MOD - 1) / n;
    mont_dot1(f, n, ivn);
  }
};

ALGO_END_NAMESPACE

#endif
