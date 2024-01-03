#define PROBLEM "https://www.luogu.com.cn/problem/P4245"

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/other/align-alloc.hpp"
#include "../../src/math/avx2/fft-radix2-twisted-avx2.hpp"

auto fft = FFT64Radix2TwistedAVX2();

constexpr u32 B = 1 << 15, X = 4; // luogu only
// constexpr u32 B = 1 << 14, X = 8;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m, P;
  fin >> n >> m >> P;
  n++, m++;
  u32 l = std::bit_ceil(n + m - 1);
  AVec<CP64> f(l * X), g(l * X);

  for (u32 i = 0; i != n; ++i) {
    u32 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      f[i * X + j] = CP64{f64(t % B)}, t /= B;
    }
  }
  for (u32 i = 0; i != m; ++i) {
    u32 t;
    fin >> t;
    for (u32 j = 0; j != X; ++j) {
      g[i * X + j] = CP64{f64(t % B)}, t /= B;
    }
  }
  FFT64Radix2TwistedAVX2 fft;
  fft.fft(f.data(), l * X);
  fft.fft(g.data(), l * X);
  fft.dot(f.data(), g.data(), l * X);
  fft.ifft(f.data(), l * X);
  fft.div2n(f.data(), l * X);
  for (u32 i = 0; i != n + m - 1; ++i) {
    u64 t = 0;
    for (u32 j = X; j != 0; --j) {
      t = (t * B + std::llround(f[i * X + j - 1].x)) % P;
    }
    fout << t << ' ';
  }
  return 0;
}