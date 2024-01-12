// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/number/mont32-const.hpp"
#include "../../src/other/align-alloc.hpp"
#include "../../src/math/avx2/ntt-radix2-twisted-avx2.hpp"
#include "../../src/number/avx2/mont-trans-helper-avx2.hpp"

constexpr std::array M = {167772161, 469762049, 754974721};
constexpr std::array G = {3, 3, 11};

template <u32 P, u32 G>
auto mul3(const AVec<u32> &f, const AVec<u32> &g, u32 l) {
  using ModT = M32C<P>;

  AVec<u32> a(l), b(l);
  auto ntt = NTT32Radix2TwistedAVX2<ModT>(G);

  auto *ma = reinterpret_cast<ModT *>(a.data());
  auto *mb = reinterpret_cast<ModT *>(b.data());

  mont32_trans_avx2(ma, f.data(), l);
  mont32_trans_avx2(mb, g.data(), l);

  ntt.ntt(ma, l);
  ntt.ntt(mb, l);
  ntt.dot(ma, mb, l);
  ntt.intt(ma, l);
  ntt.rescale(ma, l);

  mont32_get_avx2(a.data(), ma, l);
  return a;
}

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  u32 l = std::bit_ceil(n + m - 1);
  AVec<u32> f(l), g(l);
  for (u32 i = 0; i != n; ++i) {
    fin >> f[i];
  }
  for (u32 i = 0; i != m; ++i) {
    fin >> g[i];
  }

  u32 P = 1E9 + 7;
  using M32C1 = M32C<M[1]>;
  using M32C2 = M32C<M[2]>;
  u64 M12 = u64(M[0]) * M[1] % P;
  u64 inv_1 = M32C1(M[0]).inv().get();
  u64 inv_12 = (M32C2(M[0]) * M32C2(M[1])).inv().get();

  AVec<AVec<u32>> a(3);
  a[0] = mul3<M[0], G[0]>(f, g, l);
  a[1] = mul3<M[1], G[1]>(f, g, l);
  a[2] = mul3<M[2], G[2]>(f, g, l);

  for (u32 i = 0; i != n + m - 1; ++i) {
    u64 x = (a[1][i] - a[0][i] + M[1]) * inv_1 % M[1] * M[0] + a[0][i];
    u64 ans = ((a[2][i] - x % M[2] + M[2]) * inv_12 % M[2] * M12 + x) % P;
    fout << ans << ' ';
  }
  return 0;
}
