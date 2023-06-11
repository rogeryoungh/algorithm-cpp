#ifndef ALGO_MATH_POLY_INVSQRT_12ENT
#define ALGO_MATH_POLY_INVSQRT_12ENT

#include "../../base.hpp"
#include "ntt.hpp"

#include <algorithm>
#include <vector>

template <class ModT>
auto poly_invsqrt_12E(std::span<const ModT> self, u32 m, const ModT &x0) {
  u32 n = std::bit_ceil(m);
  std::vector<ModT> x(n * 2);
  x[0] = x0;
  ModT ivn2 = -ModT(2).inv();
  for (u32 t = 1; t < n; t *= 2) {
    std::vector<ModT> u(t * 4), s(t * 4);
    std::copy(self.begin(), std::min(self.begin() + t * 2, self.end()), u.begin());
    std::copy(x.begin(), x.begin() + t, s.begin());
    ntt<ModT>(u); // 4E
    ntt<ModT>(s); // 4E
    for (u32 i = 0; i < t * 4; ++i) {
      u[i] = (u[i] * s[i] * s[i] * s[i]) * ivn2;
    }
    intt<ModT>(u); // 4E
    std::copy_n(u.begin() + t, t, x.begin() + t);
  }
  return x.resize(m), x;
}

#endif // ALGO_MATH_POLY_INVSQRT_12ENT
