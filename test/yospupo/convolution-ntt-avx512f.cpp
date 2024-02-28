#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/int-only-writer.hpp"

#include "../../include/algo/math/ntt/ntt32-original-radix2-avx512f.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer<' '> fout(stdout);

  u32 n, m;
  fin >> n >> m;
  const u32 L = std::bit_ceil(n + m - 1);
  const auto M = Mont32{998244353};
  NTT32OriginalRadix2AVX512F ntt(M, 3);

  std::vector<u32> f(L), g(L);
  for (u32 i = 0; i != n; ++i)
    fin >> f[i], f[i] = M.trans(f[i]);
  for (u32 i = 0; i != m; ++i)
    fin >> g[i], g[i] = M.trans(g[i]);
  ntt.conv(f.data(), g.data(), L);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << M.get(f[i]);

  return 0;
}
