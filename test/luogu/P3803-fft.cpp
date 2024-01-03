#define PROBLEM "https://www.luogu.com.cn/problem/P3803"

#define ALGO_IO_NUMBER_ONLY

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/math/avx2/fft-radix2-twisted-avx2.hpp"
#include "../../src/other/align-alloc.hpp"

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, m;
  fin >> n >> m;
  n++, m++;
  u32 l = std::bit_ceil(n + m - 1);
  AVec<CP64> f(l);
  for (u32 i = 0; i != n; ++i) {
    u32 t;
    fin >> t, f[i].x = t;
  }
  for (u32 i = 0; i != m; ++i) {
    u32 t;
    fin >> t, f[i].y = t;
  }
  FFT64Radix2TwistedAVX2 fft;
  fft.fft(f.data(), l);
  fft.dot(f.data(), f.data(), l);
  fft.ifft(f.data(), l);
  fft.dot2(f.data(), l);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << u32(f[i].y / 2 + 0.5) << ' ';
  return 0;
}
