// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_IO_NUMBER_ONLY

#define ALGO_DISABLE_SIMD_AVX2

#include "../../src/other/fastio.hpp"
#include "../../src/math/poly/poly-base.hpp"
#include "../../src/other/modint/static-modint.hpp"

using ModT = BasicStaticModint<u32, 998244353>;
using FPS = Poly<ModT>;

i32 main() {
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
