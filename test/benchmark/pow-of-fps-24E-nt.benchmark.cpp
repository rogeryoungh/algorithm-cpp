#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#include "../../src/math/poly/exp-14E-nt-block.hpp"
#include "../../src/math/poly/pow.hpp"
#include "../../src/math/poly/ln.hpp"
#include "../../src/math/poly/div-10E-nt-block.hpp"

#include "../../src/other/modint/montgomery-space.hpp"
#include "../../src/other/modint/static-modint.hpp"
using Space = MontgomerySpace<u32, 998244353>;
using ModT = StaticModint<Space>;

/////////////////////////////////////////////////////////

#include "../benchmark-snippet.hpp"
#include <random>
struct BM_POW : TEST_BASE {
  static constexpr auto m_div = poly_div_10E_block<ModT>;
  static constexpr auto m_ln = poly_ln<ModT, m_div>;
  static constexpr auto m_exp = poly_exp_14E_block<ModT>;
  static constexpr auto m_pow = poly_pow<ModT, m_ln, m_exp>;

  void init(int n) {
    std::mt19937 rng(58);
    f.resize(n);
    f[0] = 1;
    for (u32 i = 1; i < n; i++) {
      f[i] = rng() % ModT::mod();
    }
  }
  int run(int n) {
    std::vector tf = f;
    detail::ntt_size = 0;
    m_pow(f, 998244353998244353ull, n);
    benchmark::DoNotOptimize(f[0]);
    return detail::ntt_size;
  }
  AVec<ModT> f;
};

BM_DEF(BM_POW)->RangeMultiplier(2)->Arg(1 << 18)->Arg(1 << 19)->Arg(1 << 20)->Arg(1E5)->MinTime(3);
