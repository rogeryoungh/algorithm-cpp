#ifndef ALGO_MATH_POLY_SAFE_POW
#define ALGO_MATH_POLY_SAFE_POW

#include "poly-def.hpp"
#include "../constant/prepare-inv.hpp"

template <class ModT, auto poly_pow>
auto poly_safe_pow(std::span<const ModT> f, u64 k, u64 k2, u32 m) {
  auto it = f.begin();
  while (it != f.end() && *it == 0)
    ++it;
  u32 len = it - f.begin();
  if (it == f.end() || (k > 1E9 && len >= 1) || k * len >= m) {
    AVec<ModT> r(m);
    r[0] = k == 0;
    return r;
  }
  AVec<ModT> x(it, f.end());
  ModT f0 = x[0], f0iv = f0.inv(), f0k = f0.pow(k2);
  dot_v<ModT>(x, f0iv);
  x = poly_pow(std::move(x), k, m - k * len);
  dot_v<ModT>(x, f0k);
  x.insert(x.begin(), k * len, 0);
  return x;
}

#endif // ALGO_MATH_POLY_SAFE_POW
