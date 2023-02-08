#define PROBLEM "https://judge.yosupo.jp/problem/many_aplusb"

#define ALGO_IO_NUMBER_ONLY

#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 T;
  fin >> T;
  while (T--) {
    u64 a, b;
    fin >> a >> b;
    fout << (a + b) << '\n';
  }
  return 0;
}
