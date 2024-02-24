#pragma once

#include "../base.hpp"

#include <array>
#include <cassert>
#include <concepts>
#include <cstdio>
#include <vector>

ALGO_BEGIN_NAMESPACE

template <u8 endc = 0>
struct Writer {
  std::FILE *const f;
  std::vector<u8> buf;
  std::vector<u32> pre;
  u8 *p, *end;
  Writer(std::FILE *const f, usize size = 1 << 18) : f(f), buf(size), pre(10000) {
    assert(size >= 0x100);
    p = buf.data(), end = p + size;
    for (u32 i = 0; i < 10000; ++i) {
      u32 ti = i;
      for (u32 j = 0; j != 4; ++j) {
        pre[i] = pre[i] << 8 | ti % 10 | 0x30;
        ti /= 10;
      }
    }
  }
  ~Writer() {
    flush();
  }
  void flush() {
    std::fwrite(buf.data(), 1, p - buf.data(), f);
    p = buf.data();
  }
  void reserve(usize n) {
    if (end - p < i64(n))
      flush();
  }
  void _put_u(u64 x) {
    std::array<u8, 32> tmp;
    u8 *s0 = tmp.data() + 30, *s1 = s0;
    if constexpr (endc != 0)
      *s1++ = endc;
    while (x >= 10000) {
      *reinterpret_cast<u32 *>(s0 -= 4) = pre[x % 10000];
      x /= 10000;
    }
    *reinterpret_cast<u32 *>(s0 -= 4) = pre[x % 10000];
    s0 += x < 100 ? (x < 10 ? 3 : 2) : (x < 1000 ? 1 : 0);
    p = std::copy(s0, s1, p);
  }
  template <std::integral T>
  Writer &operator<<(T x) {
    reserve(0x40);
    if (std::signed_integral<T> && x < 0) {
      *p++ = '-';
      _put_u(-x);
    } else {
      _put_u(x);
    }
    return *this;
  }
  Writer &operator<<(char x) {
    *p++ = x;
    return *this;
  }
};

ALGO_END_NAMESPACE
