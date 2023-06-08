#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_mod"

#define ALGO_IO_NUMBER_ONLY

#include "../../src/other/fastio.hpp"
#include "../../src/other/modint/dynamic-modint.hpp"
#include "../../src/other/modint/dynamic-montgomery-space.hpp"

using Space = DynamicMontgomerySpace<u32, 1>;
using ModT = DynamicModint<Space>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 T;
  fin >> T;
  while (T--) {
    u32 n, p;
    fin >> n >> p;
    ModT::set_mod(p);

    auto ret = ModT(n).sqrt();
    if (ret.has_value()) {
      fout << ret.value() << '\n';
    } else {
      fout << "-1\n";
    }
  }
  return 0;
}
