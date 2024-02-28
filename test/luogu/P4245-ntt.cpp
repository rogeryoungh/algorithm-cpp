#define PROBLEM "https://www.luogu.com.cn/problem/P4245"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/int-only-writer.hpp"

#include "../../include/algo/math/ntt/ntt32-original-radix2-avx2.hpp"
#include "../../include/algo/math/crt3-convolution.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer<' '> fout(stdout);

  u32 n, m, P;
  fin >> n >> m >> P;
  n++, m++;
  CRT3Convolution<NTT32OriginalRadix2AVX2> crt3;

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
