#ifndef ALGO_MATH_POLY_POW
#define ALGO_MATH_POLY_POW

#include "../../base.hpp"
#include "../constant/prepare-inv.hpp"

#include <span>
#include <vector>

template <static_modint_concept ModT, auto poly_ln, auto poly_exp>
auto poly_pow(std::span<const ModT> f, u64 k, u32 m) {
  ModT mk = ModT::safe(k);
  auto x = poly_ln(f, m);
  for (auto &i : x)
    i *= mk;
  return poly_exp(std::move(x), m);
}

#endif // ALGO_MATH_POLY_LN
