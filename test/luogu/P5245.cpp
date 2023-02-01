#define PROBLEM "https://www.luogu.com.cn/problem/P5245"

#include "../../src/other/fastio.hpp"
#include "../../src/math/poly/poly-base.hpp"
#include "../../src/other/modint/montgomery-space.hpp"
#include "../../src/other/modint/static-modint.hpp"

constexpr u32 P = 998244353;
using Space = MontgomerySpace<u32, P>;
using ModT = StaticModint<Space>;
using FPS = Poly<ModT>;

u32 mo(const std::string &s, u32 m) {
  u64 r = 0;
  for (const auto &si : s)
    r = (r * 10 + si - '0') % m;
  return r;
}

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n;
  std::string k;
  fin >> n >> k;
  FPS f(n);

  auto k1 = mo(k, P), k2 = mo(k, P - 1);

  for (auto &i : f)
    fin >> i;

  if (f[0] == 0 && k.size() >= 9)
    std::fill(f.begin(), f.end(), 0);
  else
    f = f.safe_pow(k1, k2, n);
  for (auto i : f)
    fout << i << ' ';
  std::cerr << std::endl << detail::ntt_size << std::endl;
  return 0;
}
