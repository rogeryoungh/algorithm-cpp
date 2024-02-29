#define PROBLEM "https://www.luogu.com.cn/problem/P3803"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/int-only-writer.hpp"
#include "../../include/algo/other/align-alloc.hpp"

#include "../../include/algo/math/ntt/ntt-original-radix2.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer<' '> fout(stdout);

  u32 n, m;
  fin >> n >> m;
  n++, m++;
  const u32 L = std::bit_ceil(n + m - 1);
  const auto M = Mont32{998244353};
  NTTOriginalRadix2 ntt(M, 3);

  AVec<u32> f(L), g(L);
  for (u32 i = 0; i != n; ++i)
    fin >> f[i];
  for (u32 i = 0; i != m; ++i)
    fin >> g[i];
  ntt.conv(f.data(), g.data(), L);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << f[i];

  return 0;
}
