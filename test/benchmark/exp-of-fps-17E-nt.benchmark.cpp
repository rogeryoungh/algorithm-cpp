#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#include "../../src/math/poly/poly-base.hpp"
#include "../../src/other/modint/basic-modint.hpp"

constexpr u32 P = 998244353;

using ModT = BasicModint<P>;
using FPS = Poly<ModT>;

/////////////////////////////////////////////////////////

#include "../benchmark-snippet.hpp"
#include <random>

struct BM_EXP : TEST_BASE {
  void init(int n) {
    std::mt19937 rng(58);
    f.resize(n);
    f[0] = 1;
    for (u32 i = 1; i < n; i++) {
      f[i] = rng() % P;
    }
  }
  int run(int n) {
    FPS tf = f;
    detail::ntt_size = 0;
    benchmark::DoNotOptimize(tf.exp(n));
    return detail::ntt_size;
  }
  FPS f;
};

BM_DEF(BM_EXP)->RangeMultiplier(2)->Arg(1 << 18)->Arg(1 << 19)->Arg(1 << 20)->Arg(1E5)->Repetitions(5);
