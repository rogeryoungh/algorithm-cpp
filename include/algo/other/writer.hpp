#pragma once

#include "../base.hpp"
#include "./itoa-helper.hpp"

#include <cassert>
#include <string>
#include <vector>

ALGO_BEGIN_NAMESPACE

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
  void _put_s(const char *s, usize n) {
    if (n >= buf.size() / 4) {
      flush();
      std::fwrite(s, 1, n, f);
    } else {
      reserve(n);
      p = std::copy(s, s + n, p);
    }
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
