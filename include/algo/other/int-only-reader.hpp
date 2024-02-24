#pragma once

#include "../base.hpp"

#include <concepts>
#include <cstdio>
#include <vector>

ALGO_BEGIN_NAMESPACE

template <class Buf>
struct Reader {
  Buf buf;
  std::vector<u16> pre;
  Reader(std::FILE *f, usize size = 1 << 18) : buf(f, size), pre(0x10000, -1) {
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
  u64 _get_u(u8 c) {
    u64 x = c & 0xf;
    while (true) {
      u16 t = *reinterpret_cast<const u16 *>(buf.p);
      auto ft = pre[t];
      buf.p += 2;
      if (ft < 100) { // len = 2
        x = x * 100 + ft;
      } else { // len = 1
        if (ft < 0x1000)
          x = x * 10 + ft - 0x100;
        else
          --buf.p;
        break;
      }
    }
    return x;
  }
  template <std::integral T>
  Reader &operator>>(T &x) {
    while (true) {
      u8 c = buf.pop();
      buf.reserve(0x40);
      if (std::signed_integral<T> && c == '-') {
        x = -T(_get_u(0));
        break;
      }
      if ('0' <= c && c <= '9') {
        x = _get_u(c);
        break;
      }
    }
    return *this;
  }
};

ALGO_END_NAMESPACE
