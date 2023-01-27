#ifndef ALGO_MATH_POLY_NT_BLOCK_HELPER
#define ALGO_MATH_POLY_NT_BLOCK_HELPER

#include "../../base.hpp"

#include <bit>
#include <span>
#include <vector>

enum { NT_BLOCK_B = 16 };

namespace detail {

std::pair<u32, u32> nt_block_len(u32 m) {
  u32 t = std::bit_ceil((m - 1) / NT_BLOCK_B + 1);
  u32 k = (m - 1) / t + 1;
  return {t, k};
}

template <class ModT>
auto nt_block_split(std::vector<ModT> &v, u32 m) {
  u32 k = v.size() / m;
  std::vector<std::span<ModT>> sv(k);
  for (u32 i = 0; i < k; ++i) {
    sv[i] = std::span<ModT>(v.begin() + i * m, m);
  }
  return sv;
}

} // namespace detail

#endif // ALGO_MATH_POLY_NT_BLOCK_HELPER
