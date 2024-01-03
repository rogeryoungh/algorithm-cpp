#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_mod"

#define ALGO_IO_NUMBER_ONLY
#define ALGO_NO_NAMESPACE

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont32-dynamic.hpp"
#include "../../src/modular/mont32-const.hpp"
#include "../../src/math/cipolla.hpp"

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 T;
  fin >> T;
  using M1 = M32D<1>;
  using M2 = M32C<998244353>;
  while (T--) {
    u32 n, p;
    fin >> n >> p;
    if (p == 998244353)
      fout << cipolla(M2(n)) << '\n';
    else {
      M1::set_mod(p);
      fout << cipolla(M1(n)) << '\n';
    }
  }
  return 0;
}
