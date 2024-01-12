// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/number/mont32-const.hpp"
#include "../../src/other/align-alloc.hpp"
#include "../../src/math/avx2/ntt-radix2-avx2.hpp"

using ModT = M32C<998244353>;
auto ntt = NTT32Radix2AVX2<ModT>(3);

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
  ntt.ntt(f.data(), l);
  ntt.ntt(g.data(), l);
  ntt.dot(f.data(), g.data(), l);
  ntt.intt(f.data(), l);
  ntt.rescale(f.data(), l);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << f[i].get() << ' ';
  return 0;
}
