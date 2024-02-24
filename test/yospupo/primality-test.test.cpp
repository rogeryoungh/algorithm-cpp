#define PROBLEM "https://judge.yosupo.jp/problem/primality_test"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/writer.hpp"

#include "../../include/algo/math/miller-rabin.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer fout(stdout);

  u32 Q;
  fin >> Q;
  while (Q--) {
    u64 n;
    fin >> n;
    fout << (miller_rabin(n) ? "Yes\n" : "No\n");
  }
  return 0;
}
