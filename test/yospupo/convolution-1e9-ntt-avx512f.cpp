#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"

#define ALGO_NO_NAMESPACE
#define ALGO_AVEC_ALIGN 64

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/int-only-writer.hpp"

#include "../../include/algo/math/ntt/ntt32-original-radix2-avx512f.hpp"
#include "../../include/algo/math/crt3-convolution.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer<' '> fout(stdout);

  u32 n, m;
  constexpr u32 P = 1E9 + 7;
  fin >> n >> m;
  CRT3Convolution<NTT32OriginalRadix2AVX512F> crt3;

  std::vector<u32> f(n), g(m);
  for (u32 i = 0; i != n; ++i)
    fin >> f[i];
  for (u32 i = 0; i != m; ++i)
    fin >> g[i];
  auto ret = crt3.conv(f, g, P);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << ret[i];

  return 0;
}
