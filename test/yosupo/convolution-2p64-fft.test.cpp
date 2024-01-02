// magic!
#pragma GCC optimize("O3")

#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_2_64"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/math/fft-radix2-twisted.hpp"

using Fft = FftR2T;
constexpr u32 B = 1 << 14, X = 16; // B * B * N <= M

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  u32 l = std::bit_ceil(n + m - 1);
  // std::vector<ModT> f(l * X), g(l * X);
  auto *f = new (std::align_val_t(32)) CP64[l * X];
  auto *g = new (std::align_val_t(32)) CP64[l * X];

  for (u32 i = 0; i != n; ++i) {
    u64 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      f[i * X + j] = CP64{f64(t % B)}, t /= B;
    }
  }
  for (u32 i = 0; i != m; ++i) {
    u64 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      g[i * X + j] = CP64{f64(t % B)}, t /= B;
    }
  }
  FftR2T fft;
  fft.fft(f, l * X);
  fft.fft(g, l * X);
  fft.dot(f, g, l * X);
  fft.ifft(f, l * X);
  fft.dot2(f, l * X);
  for (u32 i = 0; i != n + m - 1; ++i) {
    u64 t = 0;
    for (u32 j = X; j != 0; --j) {
      t = (t * B + std::llround(f[i * X + j - 1].x));
    }
    fout << t << ' ';
  }
  operator delete[](f, std::align_val_t(32));
  operator delete[](g, std::align_val_t(32));
  return 0;
}
