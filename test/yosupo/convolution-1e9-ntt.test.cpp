// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont64-const.hpp"
#include "../../src/other/align-alloc.hpp"
#include "../../src/math/ntt-radix2-twisted.hpp"

using ModT = M64C<2053641430080946177>;
NTTRadix2Twisted<ModT> ntt(7);

constexpr u32 B = 1 << 18, P = 1E9 + 7, X = 4; // B * B * N <= M

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  u32 l = std::bit_ceil(n + m - 1);

  AVec<ModT> f(l * X), g(l * X);
  for (u32 i = 0; i != n; ++i) {
    u32 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      f[i * X + j] = t % B, t /= B;
    }
  }
  for (u32 i = 0; i != m; ++i) {
    u32 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      g[i * X + j] = t % B, t /= B;
    }
  }
  ntt.ntt(f.data(), l * X);
  ntt.ntt(g.data(), l * X);
  dot(f.data(), g.data(), l * X);
  ntt.intt(f.data(), l * X);
  div2n(f.data(), l * X);
  for (u32 i = 0; i != n + m - 1; ++i) {
    u64 t = 0;
    for (u32 j = X; j != 0; --j) {
      t = (t * B + f[i * X + j - 1].get()) % P;
    }
    fout << t << ' ';
  }
  return 0;
}
