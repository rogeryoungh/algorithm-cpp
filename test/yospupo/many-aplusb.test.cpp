#define PROBLEM "https://judge.yosupo.jp/problem/many_aplusb"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/int-only-writer.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer<'\n'> fout(stdout);
  u32 T;
  fin >> T;
  while (T--) {
    u64 a, b;
    fin >> a >> b;
    fout << (a + b);
  }
  return 0;
}
