#pragma once

#include "../base.hpp"

#include <vector>

ALGO_BEGIN_NAMESPACE

struct AtoiHelper {
  std::vector<u16> pre;
  AtoiHelper() : pre(0x10000, -1) {
    for (u32 i = 0; i != 0x100; ++i) {
      for (u32 j = 0; j != 10; ++j) {
        u32 t = i * 0x100 | j | 0x30;
        if ('0' <= i && i <= '9')
          pre[t] = j * 10 + i - 0x30;
        else
          pre[t] = j | 0x100;
      }
    }
  }
  u64 getu(u8 c, const u8 *&p0) {
    const u8 *p = p0;
    u64 x = c & 0xf;
    while (true) {
      u16 t = *reinterpret_cast<const u16 *>(p);
      auto ft = pre[t];
      p += 2;
      if (ft < 100) { // len = 2
        x = x * 100 + ft;
      } else { // len = 1
        if (ft < 0x1000)
          x = x * 10 + ft - 0x100;
        else
          --p;
        break;
      }
    }
    return p0 = p, x;
  }
};

ALGO_END_NAMESPACE
