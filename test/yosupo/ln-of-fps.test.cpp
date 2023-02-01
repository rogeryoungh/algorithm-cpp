// magic!
#pragma GCC target("abm,movbe,bmi,bmi2,lzcnt,popcnt,avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("inline")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("omit-frame-pointer")

#define PROBLEM "https://judge.yosupo.jp/problem/log_of_formal_power_series"

#include "../../src/other/fastio.hpp"
#include "../../src/other/modint/static-modint.hpp"
#include "../../src/math/poly/poly-base.hpp"

using ModT = BasicStaticModint<u32, 998244353>;
using FPS = Poly<ModT>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n;
  fin >> n;
  FPS f(n);

  for (auto &i : f)
    fin >> i;
  for (auto i : f.ln(n))
    fout << i << ' ';
  std::cerr << std::endl << detail::ntt_size << std::endl;
  return 0;
}
