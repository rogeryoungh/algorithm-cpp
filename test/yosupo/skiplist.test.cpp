#define PROBLEM "https://judge.yosupo.jp/problem/associative_array"

#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/datastruct/skiplist.hpp"

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  // auto &fin = std::cin;
  // auto &fout = std::cout;
  i32 Q;
  fin >> Q;
  SkipList<u64, u64> skiplist;
  while (Q--) {
    u32 o;
    fin >> o;
    if (o == 0) {
      u64 k, v;
      fin >> k >> v;
      skiplist[k] = v;
    } else {
      u64 k;
      fin >> k;
      fout << skiplist[k] << '\n';
    }
  }
  return 0;
}