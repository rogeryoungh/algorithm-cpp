// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_IO_NUMBER_ONLY

#define ALGO_DYNAMIC_MOD

#include "../../src/other/fastio.hpp"
#include "../../src/math/poly/poly-base.hpp"

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

i32 main() {
#ifdef ALGO_DYNAMIC_MOD
  ModT::set_mod(998244353);
#endif
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  FPS f(n), g(m);

  for (auto &i : f)
    fin >> i;
  for (auto &i : g)
    fin >> i;
  for (auto i : f *g)
    fout << i << ' ';
  std::cerr << std::endl << detail::ntt_size << std::endl;
  return 0;
}
