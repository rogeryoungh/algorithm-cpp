#define PROBLEM "https://loj.ac/p/143"

#include "../../src/math/primality-test/baillie-psw.hpp"
#include "../../src/other/fastio.hpp"

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  while (true) {
    u64 n;
    fin >> n;
    if (n == 0)
      break;
    fout << (baillie_psw_test(n) ? 'Y' : 'N') << '\n';
  }
  return 0;
}
