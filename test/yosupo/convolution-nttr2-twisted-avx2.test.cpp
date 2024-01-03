// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont32-const.hpp"
#include "../../src/math/avx2/ntt-radix2-twisted-avx2.hpp"

using ModT = M32C<998244353>;
using NTT = NTT32Radix2TwistedAVX2<ModT, 3>;

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
  NTT::set_mod();
  NTT::ntt(f.data(), l);
  NTT::ntt(g.data(), l);
  NTT::dot(f.data(), g.data(), l);
  NTT::intt(f.data(), l);
  NTT::dot2(f.data(), l);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << f[i].get() << ' ';
  return 0;
}
