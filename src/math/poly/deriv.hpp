#ifndef ALGO_MATH_POLY_DERIV
#define ALGO_MATH_POLY_DERIV

#include "poly-def.hpp"
#include "../constant/prepare-inc.hpp"

template <class ModT>
auto poly_deriv(std::span<const ModT> f, u32 m) {
  const auto &inc = prepare_inc<ModT>(m);
  AVec<ModT> x(m);
  std::copy(f.begin() + 1, std::min(f.begin() + m + 1, f.end()), x.begin());
  dot<ModT>(x, {inc.begin() + 1, m});
  return x;
}

#endif // ALGO_MATH_POLY_DERIV
