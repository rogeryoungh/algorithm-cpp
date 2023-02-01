#define PROBLEM "https://loj.ac/p/143"

#include "../../src/other/fastio.hpp"
#include "../../src/other/modint/static-modint.hpp"
#include "../../src/math/poly/poly-base.hpp"
#include "../../src/other/modint/montgomery-space.hpp"

using Space = MontgomerySpace<u32, 998244353>;
using ModT = StaticModint<Space>;
using FPS = Poly<ModT>;

i32 main() {
  FastI fin(stdin);
  FastO fout(stdout);
  u32 n, k;
  fin >> n >> k;
  n++;
  FPS f(n);
  for (auto &i : f)
    fin >> i;
  auto ans = ((f - FPS{f[0] - 2} - f.sqrt(n).inv(n).integr(n).exp(n)).ln(n) + FPS{ModT(1)}).pow(k, n).deriv(n - 1);
  for (const auto &i : ans)
    fout << i << ' ';
  return 0;
}
