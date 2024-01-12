#define PROBLEM "https://www.luogu.com.cn/problem/P3803"

#define ALGO_IO_NUMBER_ONLY

#define ALGO_NO_NAMESPACE
#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/number/mont32-dynamic.hpp"
#include "../../src/other/align-alloc.hpp"
#include "../../src/math/avx2/ntt-radix2-twisted-avx2.hpp"
#include "../../src/number/avx2/mont-io-helper-avx2.hpp"

using ModT = M32D<1>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  ModT::set_mod(998244353);
  auto ntt = NTT32Radix2TwistedAVX2<ModT>(3);

  u32 n, m;
  fin >> n >> m;
  n++, m++;
  u32 l = std::bit_ceil(n + m - 1);

  AVec<ModT> f(l), g(l);
  mont_read(fin, f.data(), n);
  mont_read(fin, g.data(), m);
  ntt.ntt(f.data(), l);
  ntt.ntt(g.data(), l);
  ntt.dot(f.data(), g.data(), l);
  ntt.intt(f.data(), l);
  ntt.rescale(f.data(), l);
  mont_write(fout, f.data(), n + m - 1);
  return 0;
}
