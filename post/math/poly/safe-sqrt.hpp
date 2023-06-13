#ifndef ALGO_MATH_POLY_SAFE_SQRT
#define ALGO_MATH_POLY_SAFE_SQRT

#include "poly-def.hpp"
#include "../../other/align-alloc.hpp"

#include <optional>

template <class ModT, auto poly_sqrt>
std::optional<AVec<ModT>> poly_safe_sqrt(std::span<const ModT> f, u32 m) {
  auto it = f.begin();
  while (it != f.end() && *it == 0)
    ++it;
  if (it == f.end())
    return AVec<ModT>(m);
  auto sq = it->sqrt();
  u32 len = it - f.begin();
  if (len % 2 == 1 || !sq.has_value())
    return std::nullopt;
  AVec<ModT> x = {it, f.end()};
  x = poly_sqrt(x, m - len / 2, sq.value());
  x.insert(x.begin(), len / 2, 0);
  return x;
}

#endif // ALGO_MATH_POLY_SAFE_SQRT
