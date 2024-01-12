#ifndef ALGO_H_MATH_MONT32_IO_HELPER
#define ALGO_H_MATH_MONT32_IO_HELPER

#include "../base.hpp"

ALGO_BEGIN_NAMESPACE

template <class ModT>
void mont_read(auto &io, ModT *f, u32 n) {
  for (u32 i = 0; i != n; ++i) {
    u32 t;
    io >> t;
    f[i] = t;
  }
}

template <class ModT>
void mont_write(auto &io, ModT *f, u32 n) {
  io << f[0].get();
  for (u32 i = 1; i != n; ++i) {
    io << ' ' << f[i].get();
  }
}

ALGO_END_NAMESPACE

#endif
