// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/other/align-alloc.hpp"
#include "../../src/math/fft-radix2-twisted.hpp"

FFTRadix2Twisted fft;

constexpr u32 B = 1 << 14, P = 1E9 + 7, X = 8; // B * B * N <= M

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  u32 l = std::bit_ceil(n + m - 1);

  AVec<CP64> f(l * X);
  for (u32 i = 0; i != n; ++i) {
    u32 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      f[i * X + j].x = t % B, t /= B;
    }
  }
  for (u32 i = 0; i != m; ++i) {
    u32 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      f[i * X + j].y = t % B, t /= B;
    }
  }
  fft.fft(f.data(), l * X);
  fft.dot(f.data(), f.data(), l * X);
  fft.ifft(f.data(), l * X);
  fft.rescale(f.data(), l * X);
  for (u32 i = 0; i != n + m - 1; ++i) {
    u64 t = 0;
    for (u32 j = X; j != 0; --j) {
      t = (t * B + std::llround(f[i * X + j - 1].y / 2)) % P;
    }
    fout << t << ' ';
  }
  return 0;
}
