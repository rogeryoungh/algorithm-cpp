#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

// #define ALGO_DISABLE_SIMD_AVX2
// #define ALGO_DISABLE_NTT_RADIX_4
// #define ALGO_DISABLE_NTT_CLASSICAL

// #define ALGO_DYNAMIC_MOD

#include "../../src/math/poly/poly-base.hpp"

constexpr u32 P = 998244353;

#ifndef ALGO_DYNAMIC_MOD
#include "../../src/other/modint/montgomery-space.hpp"
#include "../../src/other/modint/static-modint.hpp"
using Space = MontgomerySpace<u32, 998244353>;
using ModT = StaticModint<Space>;
#else
#include "../../src/other/modint/dynamic-montgomery-space.hpp"
#include "../../src/other/modint/dynamic-modint.hpp"
using Space = DynamicMontgomerySpace<u32, 1>;
using ModT = DynamicModint<Space>;
#endif

using FPS = Poly<ModT>;

/////////////////////////////////////////////////////////

#include "../benchmark-snippet.hpp"
#include <random>

struct BM_MUL : TEST_BASE {
  void init(int n) {
#ifdef ALGO_DYNAMIC_MOD
    ModT::set_mod(P);
#endif
    std::mt19937 rng(58);
    f.resize(n), g.resize(n);
    for (u32 i = 0; i < n; i++) {
      f[i] = rng() % P;
      g[i] = rng() % P;
    }
  }
  int run(int) {
    FPS tf = f, tg = g;
    detail::ntt_size = 0;
    benchmark::DoNotOptimize(tf * tg);
    return detail::ntt_size;
  }
  FPS f, g;
};

BM_DEF(BM_MUL)->RangeMultiplier(2)->Arg(1 << 18)->Arg(1 << 19)->Arg(1 << 20)->Arg(1E5)->MinTime(3);
