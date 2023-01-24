#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_formal_power_series"

#include "../../src/other/fastio.hpp"
#include "../../src/other/modint/basic-modint.hpp"
#include "../../src/math/poly/poly-base.hpp"

using ModT = BasicModint<998244353>;
using FPS = Poly<ModT>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n;
  fin >> n;
  FPS f(n);

  for (auto &i : f)
    fin >> i;
  for (auto i : f.inv(n))
    fout << i << ' ';
  std::cerr << std::endl << detail::ntt_size << std::endl;
  return 0;
}
