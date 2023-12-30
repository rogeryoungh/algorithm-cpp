// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont32-const.hpp"
#include "../../src/math/avx512f/ntt-radix2-avx512f.hpp"

using ModT = M32C<998244353>;
using Ntt = Ntt32R2Avx512F<ModT, 3>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  u32 l = std::bit_ceil(n + m - 1);
  auto *f = new (std::align_val_t(64)) ModT[l];
  auto *g = new (std::align_val_t(64)) ModT[l];
  for (u32 i = 0; i != n; ++i) {
    u32 t;
    fin >> t, f[i] = t;
  }
  for (u32 i = 0; i != m; ++i) {
    u32 t;
    fin >> t, g[i] = t;
  }
  Ntt::setMod();
  Ntt::ntt(f, l);
  Ntt::ntt(g, l);
  Ntt::dot(f, g, l);
  Ntt::intt(f, l);
  Ntt::dot2(f, l);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << f[i].get() << ' ';
  operator delete[](f, std::align_val_t(64));
  operator delete[](g, std::align_val_t(64));
  return 0;
}
