#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_mod"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/int-only-writer.hpp"

#include "../../include/algo/math/cipolla.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer<'\n'> fout(stdout);
  u32 T;
  fin >> T;
  auto M9 = Mont32{998244353};
  while (T--) {
    u32 n, p;
    fin >> n >> p;
    if (p == 998244353) {
      fout << cipolla(M9, n);
    } else {
      fout << cipolla(Mont32{p}, n);
    }
  }
  return 0;
}
