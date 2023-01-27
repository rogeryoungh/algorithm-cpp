#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#include "../../src/math/poly/poly-base.hpp"
#include "../../src/other/modint/basic-modint.hpp"

constexpr u32 P = 998244353;

using ModT = BasicModint<P>;
using FPS = Poly<ModT>;

/////////////////////////////////////////////////////////

#include "../benchmark-snippet.hpp"
#include <random>

struct BM_MUL : TEST_BASE {
  void init(int n) {
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
