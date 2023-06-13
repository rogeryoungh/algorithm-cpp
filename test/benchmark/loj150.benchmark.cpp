#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

// #define ALGO_DISABLE_SIMD_AVX2
// #define ALGO_DISABLE_NTT_RADIX_4
// #define ALGO_DISABLE_NTT_CLASSICAL

#include "../../src/math/poly/poly-base.hpp"

#include "../../src/other/modint/montgomery-space.hpp"
#include "../../src/other/modint/static-modint.hpp"
using Space = MontgomerySpace<u32, 998244353>;
using ModT = StaticModint<Space>;

using FPS = Poly<ModT>;

/////////////////////////////////////////////////////////

#include "../benchmark-snippet.hpp"
#include <random>

struct BM_SQRT : TEST_BASE {
  void init(int n) {
    std::mt19937 rng(58);
    f.resize(n);
    f[0] = 1;
    for (u32 i = 1; i < n; i++) {
      f[i] = rng() % ModT::mod();
    }
  }
  int run(int n) {
    AVec<ModT> tf = f;
    detail::ntt_size = 0;
    auto ans = ((f - FPS{f[0] - 2} - f.invsqrt(n).integr(n).exp(n)).ln(n) + FPS{ModT(1)}).pow(19260817, n).deriv(n - 1);

    benchmark::DoNotOptimize(ans);
    return detail::ntt_size;
  }
  FPS f;
};

BM_DEF(BM_SQRT)->RangeMultiplier(2)->Arg(1 << 16)->Arg(1 << 17)->Arg(1 << 18)->Arg(1E5)->MinTime(3);
