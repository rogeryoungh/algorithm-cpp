#define PROBLEM "https://www.luogu.com.cn/problem/P3803"

#include "../../src/other/fastio.hpp"
#include "../../src/other/modint/basic-modint.hpp"
#include "../../src/math/poly/poly-base.hpp"

using ModT = BasicModint<998244353>;
using FPS = Poly<ModT>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  n++, m++;
  FPS f(n), g(m);

  for (auto &i : f)
    fin >> i;
  for (auto &i : g)
    fin >> i;
  for (auto i : f *g)
    fout << i << ' ';
  std::cerr << std::endl << ntt_size << std::endl;
  return 0;
}
