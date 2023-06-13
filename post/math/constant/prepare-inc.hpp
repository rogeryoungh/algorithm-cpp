#ifndef ALGO_MATH_CONSTANT_PREPARE_INC
#define ALGO_MATH_CONSTANT_PREPARE_INC

#include "../../base.hpp"

#include "../../other/align-alloc.hpp"

template <class ModT>
auto &prepare_inc(u32 m) {
  static AVec<ModT> inc{0, 1};
  m = std::bit_ceil(m);
  if (inc.size() < m) {
    u32 n = inc.size();
    inc.resize(m);
    for (u32 i = n; i < m; ++i) {
      inc[i] = i;
    }
  }
  return inc;
}

#endif // ALGO_MATH_CONSTANT_PREPARE_INC
