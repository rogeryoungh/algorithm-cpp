#define PROBLEM "https://www.luogu.com.cn/problem/P5491"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/writer.hpp"

#include "../../include/algo/math/cipolla.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer fout(stdout);
  u32 T;
  fin >> T;
  while (T--) {
    u32 n, p;
    fin >> n >> p;
    i32 ret = cipolla(Mont32{p}, n);
    if (ret == -1) {
      fout << "Hola!\n";
    } else if (ret * 2u == p) {
      fout << ret << '\n';
    } else {
      fout << ret << ' ' << p - ret << '\n';
    }
  }
  return 0;
}
