#define PROBLEM "https://www.luogu.com.cn/problem/P3803"

#define ALGO_IO_NUMBER_ONLY
// #define ALGO_DISABLE_NTT_RADIX_4
// #define ALGO_DISABLE_SIMD_AVX2
// #define ALGO_DISABLE_NTT_CLASSICAL
// #define ALGO_DYNAMIC_MOD

#include "../../src/other/fastio.hpp"
#include "../../src/math/poly/ntt.hpp"
#include "../../src/math/poly/vec-dots.hpp"
#include "../../src/other/align-alloc.hpp"

#ifndef ALGO_DYNAMIC_MOD
#include "../../src/other/modint/montgomery-space.hpp"
#include "../../src/other/modint/static-modint.hpp"
using Space = MontgomerySpace<u32, 998244353>;
using ModT = StaticModint<Space>;
#else
#include "../../src/other/modint/dynamic-montgomery-space.hpp"
#include "../../src/other/modint/dynamic-modint.hpp"
using Space = DynamicMontgomerySpace<u32, 1>;
using ModT = DynamicModint<Space>;
#endif

i32 main() {
#ifdef ALGO_DYNAMIC_MOD
  ModT::set_mod(998244353);
#endif
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  n++, m++;
  u32 L = std::bit_ceil(n + m - 1);
  AVec<ModT> f(L), g(L);
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
