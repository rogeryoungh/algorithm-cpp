#ifndef ALGO_MATH_POLY_POW
#define ALGO_MATH_POLY_POW

#include "poly-def.hpp"
#include "../constant/prepare-inv.hpp"

template <class ModT, auto poly_ln, auto poly_exp>
auto poly_pow(std::span<const ModT> f, u64 k, u32 m) {
  if (k == 0) {
    AVec<ModT> x(m);
    x[0] = 1;
    return x;
  } else {
    ModT mk = ModT::safe(k);
    auto x = poly_ln(f, m);
    dot_v<ModT>(x, mk);
    return poly_exp(std::move(x), m);
  }
}

#endif // ALGO_MATH_POLY_LN
