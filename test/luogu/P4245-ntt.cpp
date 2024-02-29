#define PROBLEM "https://www.luogu.com.cn/problem/P4245"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/int-only-writer.hpp"

#if 1
#define ALGO_AVEC_ALIGN 64
#include "../../include/algo/math/ntt/ntt32-original-radix2-avx512f.hpp"
using NTT = NTT32OriginalRadix2AVX512F;
#else
#include "../../include/algo/math/ntt/ntt32-original-radix2-avx2.hpp"
using NTT = NTT32OriginalRadix2AVX2;
#endif

#include "../../include/algo/math/crt3-convolution.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer<' '> fout(stdout);

  u32 n, m, P;
  fin >> n >> m >> P;
  n++, m++;

  AVec<u32> f(n), g(m);
  for (u32 i = 0; i != n; ++i)
    fin >> f[i];
  for (u32 i = 0; i != m; ++i)
    fin >> g[i];
  auto ret = convolution_crt3<NTT>(f, g, P);
  for (u32 i = 0; i != n + m - 1; ++i)
    fout << ret[i];

  return 0;
}
