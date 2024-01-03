#define PROBLEM "https://www.luogu.com.cn/problem/P5491"

#define ALGO_IO_NUMBER_ONLY
#define ALGO_NO_NAMESPACE

#include "../../src/other/fastio.hpp"
#include "../../src/modular/mont32-dynamic.hpp"
#include "../../src/math/cipolla.hpp"

using ModT = M32D<1>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 T;
  fin >> T;
  while (T--) {
    u32 n, p;
    fin >> n >> p;
    ModT::set_mod(p);
    auto ret = cipolla(ModT(n));
    if (ret != -1) {
      if (ret == -ret)
        fout << ret << '\n';
      else
        fout << ret << ' ' << p - ret << '\n';
    } else {
      fout << "Hola!\n";
    }
  }
  return 0;
}
