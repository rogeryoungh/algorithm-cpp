#define PROBLEM "https://www.luogu.com.cn/problem/P3803"

#define ALGO_IO_NUMBER_ONLY

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont32-const.hpp"
#include "../../src/math/ntt-radix2-twisted.hpp"

using ModT = M32C<998244353>;
using Ntt = NttR2T<ModT, 3>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  n++, m++;
  u32 l = std::bit_ceil(n + m - 1);
  auto *f = new (std::align_val_t(32)) ModT[l];
  auto *g = new (std::align_val_t(32)) ModT[l];
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
  operator delete[](f, std::align_val_t(32));
  operator delete[](g, std::align_val_t(32));
  return 0;
}
