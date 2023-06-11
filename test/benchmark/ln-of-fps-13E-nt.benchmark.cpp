#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#include "../../src/math/poly/div-13E-nt.hpp"
#include "../../src/math/poly/ln.hpp"

#include "../../src/other/modint/montgomery-space.hpp"
#include "../../src/other/modint/static-modint.hpp"
using Space = MontgomerySpace<u32, 998244353>;
using ModT = StaticModint<Space>;

/////////////////////////////////////////////////////////

#include "../benchmark-snippet.hpp"
#include <random>

struct BM_LN : TEST_BASE {
  void init(int n) {
    std::mt19937 rng(58);
    f.resize(n);
    f[0] = 1;
    for (u32 i = 1; i < n; i++) {
      f[i] = rng() % ModT::mod();
    }
  }
  int run(int n) {
    auto tf = f;
    detail::ntt_size = 0;
    poly_ln<ModT, poly_div_13E<ModT>>(tf, n);
    benchmark::DoNotOptimize(tf[0]);
    return detail::ntt_size;
  }
  AVec<ModT> f;
};

BM_DEF(BM_LN)->RangeMultiplier(2)->Arg(1 << 18)->Arg(1 << 19)->Arg(1 << 20)->Arg(1E5)->MinTime(3);
