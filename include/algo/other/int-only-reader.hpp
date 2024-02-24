#pragma once

#include "../base.hpp"
#include "./atoi-helper.hpp"

#include <cstdio>

ALGO_BEGIN_NAMESPACE

template <class Buf>
struct Reader {
  Buf buf;
  AtoiHelper atoi;
  Reader(std::FILE *f, usize size = 1 << 18) : buf(f, size) {}
  bool eof() const {
    return buf.eof();
  }
  template <std::integral T>
  Reader &operator>>(T &x) {
    while (true) {
      buf.reserve(0x40);
      u8 c = buf.pop();
      if (std::signed_integral<T> && c == '-') {
        x = -T(atoi.getu(0, buf.p));
        break;
      }
      if ('0' <= c && c <= '9') {
        x = atoi.getu(c, buf.p);
        break;
      }
    }
    return *this;
  }
};

ALGO_END_NAMESPACE
