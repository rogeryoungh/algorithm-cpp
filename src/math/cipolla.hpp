#ifndef ALGO_H_MATH_CIPOLLA
#define ALGO_H_MATH_CIPOLLA

#include "../base.hpp"
#include <algorithm>

ALGO_BEGIN_NAMESPACE

template <class ModT>
ModT legendre(ModT a) {
  return a.pow((ModT::MOD - 1) / 2);
}

template <class ModT>
i64 cipolla(ModT n) {
  if (n == 0)
    return 0;
  if (legendre(n) != 1)
    return -1;
  if (ModT::MOD == 2)
    return 1;
  for (u64 a = 0; a < ModT::MOD; a++) {
    ModT i = a;
    i = i * i - n;
    using FP2 = std::pair<ModT, ModT>;
    auto mul = [i](const FP2 &l, const FP2 &r) {
      auto [la, lb] = l;
      auto [ra, rb] = r;
      return FP2{la * ra + lb * rb * i, lb * ra + la * rb};
    };
    if (legendre(i) == -ModT{1}) {
      FP2 x = {1, 1}, u = {a, 1};
      for (u64 b = (ModT::MOD + 1) / 2; b; b /= 2) {
        if (b % 2 == 1)
          x = mul(x, u);
        u = mul(u, u);
      }
      return std::min(x.first.get(), (-x.first).get());
    }
  }
  return -1;
}

ALGO_END_NAMESPACE

#endif
