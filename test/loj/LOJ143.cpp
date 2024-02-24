#define PROBLEM "https://loj.ac/p/143"

#define ALGO_NO_NAMESPACE

#include "../../include/algo/other/mmap-buffer.hpp"
#include "../../include/algo/other/int-only-reader.hpp"
#include "../../include/algo/other/writer.hpp"

#include "../../include/algo/math/miller-rabin.hpp"

i32 main() {
  Reader<MmapBuf> fin(stdin);
  Writer fout(stdout);

  while (!fin.eof()) {
    u64 n;
    fin >> n;
    fout << (miller_rabin(n) ? "Y\n" : "N\n");
  }
  return 0;
}
