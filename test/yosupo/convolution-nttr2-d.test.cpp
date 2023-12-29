// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont32-dynamic.hpp"
#include "../../src/math/ntt-radix2.hpp"

using ModT = M32D<1>;
using Ntt = NttR2<ModT, 3>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  ModT::setMod(998244353);
  u32 l = std::bit_ceil(n + m - 1);
  std::vector<ModT> f(l), g(l);
  for (u32 i = 0; i != n; ++i) {
    u32 t;
    fin >> t, f[i] = t;
  }
  for (u32 i = 0; i != m; ++i) {
    u32 t;
    fin >> t, g[i] = t;
  }
  Ntt::setMod();
  Ntt::ntt(f.data(), l);
  Ntt::ntt(g.data(), l);
  Ntt::dot(f.data(), g.data(), l);
  Ntt::intt(f.data(), l);
  Ntt::dot2(f.data(), l);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << f[i].get() << ' ';
  return 0;
}
