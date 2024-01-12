#ifndef ALGO_H_NUMBER_AVX2_MONT32_IO_HELPER
#define ALGO_H_NUMBER_AVX2_MONT32_IO_HELPER

#include "./mont32x8.hpp"
#include <array>

ALGO_BEGIN_NAMESPACE

template <class ModT>
void mont_read(auto &io, ModT *f, u32 n) {
  alignas(i256) std::array<u32, 8> r;
  u32 i = 0;
  for (; i + 7 < n; i += 8) {
    for (u32 j = 0; j != 8; ++j)
      io >> r[j];
    u32x8 rx = i256_load(r.data());
    _mm256_store_si256(reinterpret_cast<u32x8 *>(f + i), M32x8<ModT>::trans(rx));
  }
  for (; i != n; ++i) {
    u32 t;
    io >> t;
    f[i] = t;
  }
}

template <class ModT>
void mont_write(auto &io, ModT *f, u32 n) {
  alignas(i256) std::array<u32, 8> r;
  u32 i = 0;
  for (; i + 8 < n; i += 8) {
    u32x8 rx = i256_load(f + i);
    _mm256_store_si256(reinterpret_cast<u32x8 *>(r.data()), M32x8<ModT>::get(rx));
    for (u32 j = 0; j != 8; ++j) {
      io << r[j] << ' ';
    }
  }
  io << f[i++].get();
  for (; i != n; ++i) {
    io << ' ' << f[i].get();
  }
}

ALGO_END_NAMESPACE

#endif
