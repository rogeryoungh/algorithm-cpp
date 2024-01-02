// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont32-const.hpp"
#include "../../src/math/avx2/ntt-radix4-avx2.hpp"
#include "../../src/other/align-alloc.hpp"

using ModT = M32C<998244353>;
using Ntt = Ntt32R4Avx2<ModT, 3>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  u32 l = std::bit_ceil(n + m - 1);
  AVec<ModT> f(l), g(l);
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
