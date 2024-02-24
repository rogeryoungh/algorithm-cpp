#pragma once

#include "../base.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <string>
#include <vector>

ALGO_BEGIN_NAMESPACE

struct Writer {
  std::FILE *const f;
  std::vector<u8> buf;
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
  void _put_u(u64 x) {
    std::array<u8, 32> tmp;
    u8 *s0 = tmp.data() + 30, *s1 = s0;
    do {
      *--s0 = x % 10 + '0', x /= 10;
    } while (x > 0);
    p = std::copy(s0, s1, p);
  }
  void _put_s(const char *s, usize n) {
    if (n >= buf.size() / 4) {
      flush();
      std::fwrite(s, 1, n, f);
    } else {
      reserve(n);
      p = std::copy_n(s, n, p);
    }
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
    reserve(0x40);
    *p++ = x;
    return *this;
  }
  Writer &operator<<(const std::string &x) {
    _put_s(x.data(), x.size());
    return *this;
  }
  template <usize N>
  Writer &operator<<(const char (&s)[N]) {
    _put_s(s, N - 1);
    return *this;
  }
};

ALGO_END_NAMESPACE
