#pragma once

#include "../base.hpp"

#include <cstring>
#include <vector>
#include <array>

ALGO_BEGIN_NAMESPACE

struct ItoaHelper {
  std::vector<u32> pre;
  ItoaHelper() : pre(10000) {
    for (u32 i = 0; i < 10000; ++i) {
      u32 ti = i;
      for (u32 j = 0; j != 4; ++j) {
        pre[i] = pre[i] << 8 | ti % 10 | 0x30;
        ti /= 10;
      }
    }
  }
  void putu(u64 x, u8 *&p) {
    std::array<u8, 32> tmp;
    u8 *s0 = tmp.data() + 30, *s1 = s0;
    while (x >= 10000) {
      std::memcpy(s0 -= 4, &pre[x % 10000], 4);
      x /= 10000;
    }
    std::memcpy(s0 -= 4, &pre[x % 10000], 4);
    s0 += x < 100 ? (x < 10 ? 3 : 2) : (x < 1000 ? 1 : 0);
    p = std::copy(s0, s1, p);
  }
};

ALGO_END_NAMESPACE
