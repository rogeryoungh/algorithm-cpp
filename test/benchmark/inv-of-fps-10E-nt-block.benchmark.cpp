#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#include "../../src/math/poly/inv-10E-nt-block.hpp"

#include "../../src/other/modint/static-modint.hpp"

constexpr u32 P = 998244353;

using ModT = BasicStaticModint<u32, P>;

/////////////////////////////////////////////////////////

#include "../benchmark-snippet.hpp"
#include <random>

struct BM_INV : TEST_BASE {
  void init(int n) {
    std::mt19937 rng(58);
    f.resize(n);
    f[0] = rng() % (P - 1) + 1;
    for (u32 i = 1; i < n; i++) {
      f[i] = rng() % P;
    }
  }
  int run(int n) {
    auto tf = f;
    detail::ntt_size = 0;
    poly_inv_10E_block<ModT>(tf, n);
    benchmark::DoNotOptimize(tf[0]);
    return detail::ntt_size;
  }
  std::vector<ModT> f;
};

BM_DEF(BM_INV)->RangeMultiplier(2)->Arg(1 << 18)->Arg(1 << 19)->Arg(1 << 20)->Arg(1E5)->MinTime(3);
