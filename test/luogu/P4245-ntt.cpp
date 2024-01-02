#define PROBLEM "https://www.luogu.com.cn/problem/P4245"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont64-const.hpp"
#include "../../src/math/ntt-radix4.hpp"

using ModT = M64C<2053641430080946177>;
using Ntt = NttR4<ModT, 7>;

constexpr u32 B = 1 << 18, X = 4; // B * B * N <= M

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m, P;
  fin >> n >> m >> P;
  n++, m++;
  u32 l = std::bit_ceil(n + m - 1);
  // std::vector<ModT> f(l * X), g(l * X);
  auto *f = new (std::align_val_t(32)) ModT[l * X];
  auto *g = new (std::align_val_t(32)) ModT[l * X];

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
  Ntt::setMod();
  Ntt::ntt(f, l * X);
  Ntt::ntt(g, l * X);
  Ntt::dot(f, g, l * X);
  Ntt::intt(f, l * X);
  Ntt::dot2(f, l * X);
  for (u32 i = 0; i != n + m - 1; ++i) {
    u64 t = 0;
    for (u32 j = X; j != 0; --j) {
      t = (t * B + f[i * X + j - 1].get()) % P;
    }
    fout << t << ' ';
  }
  operator delete[](f, std::align_val_t(32));
  operator delete[](g, std::align_val_t(32));
  return 0;
}
