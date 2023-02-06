#define PROBLEM "https://www.luogu.com.cn/problem/P3803"

#include "../../src/other/fastio.hpp"
#include "../../src/other/modint/montgomery-space.hpp"
#include "../../src/other/modint/static-modint.hpp"
#include "../../src/math/poly/ntt.hpp"
#include "../../src/math/poly/vec-dots.hpp"

using Space = MontgomerySpace<u32, 998244353>;
using ModT = StaticModint<Space>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  n++, m++;
  u32 L = std::bit_ceil(n + m - 1);
  std::vector<ModT> f(L), g(L);
  for (i32 i = 0; i < n; ++i)
    fin >> f[i];
  for (i32 i = 0; i < m; ++i)
    fin >> g[i];
  ntt<ModT>(f), ntt<ModT>(g);
  dot<ModT>(f, g);
  intt<ModT>(f);
  for (i32 i = 0; i < n + m - 1; ++i)
    fout << f[i] << ' ';
  std::cerr << std::endl << detail::ntt_size << std::endl;
  return 0;
}
