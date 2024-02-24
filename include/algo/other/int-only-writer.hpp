#pragma once

#include "../base.hpp"
#include "./itoa-helper.hpp"

#include <cassert>
#include <cstdio>
#include <vector>

ALGO_BEGIN_NAMESPACE

template <u8 endc = 0>
struct Writer {
  std::FILE *const f;
  std::vector<u8> buf;
  ItoaHelper itoa;
  u8 *p, *end;
  Writer(std::FILE *const f, usize size = 1 << 18) : f(f), buf(size) {
    assert(size >= 0x100);
    p = buf.data(), end = p + size;
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
  template <std::integral T>
  Writer &operator<<(T x) {
    reserve(0x40);
    if (std::signed_integral<T> && x < 0) {
      *p++ = '-';
      itoa.putu(-x, p);
    } else {
      itoa.putu(x, p);
    }
    if constexpr (endc != 0)
      *p++ = endc;
    return *this;
  }
  Writer &operator<<(char x) {
    *p++ = x;
    return *this;
  }
};

ALGO_END_NAMESPACE
