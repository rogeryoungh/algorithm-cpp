#ifndef ALGO_FASTIO
#define ALGO_FASTIO

#include "../base.hpp"

#include <cstring>
#include <string>
#include <vector>
#include <array>

// 快读，未实现 string 读入，未限定读入类型

struct BasicBuffer {
  std::vector<char> s;
  FILE *f;
  BasicBuffer(FILE *f_, u32 sz = 1 << 18) : s(sz), f(f_) {}
  char *p = s.data(), *beg = p, *end = p + s.size();
  inline char getc() {
    if (p == end) {
      u32 r = std::fread(beg, 1, end - beg, f);
      std::memset(beg + r, 0, s.size() - r);
      p = s.data();
    }
    return *p++;
  }
  inline void putc(char c) {
    if (p == end)
      flush();
    *p++ = c;
  }
  inline void puts(const char *x) {
    while (*x != 0)
      putc(*x++);
  }
  void flush() {
    std::fwrite(beg, 1, p - beg, f);
    p = s.data();
  }
};

struct FastI : BasicBuffer {
  using BasicBuffer::BasicBuffer;
  FastI(FILE *f, u32 sz = 1 << 18) : BasicBuffer(f, sz) {
    std::fread(beg, 1, end - beg, f);
  }
  FastI &operator>>(char &x) {
    return x = getc(), self;
  }
  template <std::unsigned_integral T>
  FastI &operator>>(T &x) {
    x = 0;
    char c = getc();
    while (!std::isdigit(c) && c != 0)
      c = getc();
    while (std::isdigit(c))
      x = x * 10 + c - '0', c = getc();
    return self;
  }
  template <std::signed_integral T>
  FastI &operator>>(T &x) {
    x = 0;
    char c = getc();
    bool sgn = true;
    while (!std::isdigit(c) && c != 0)
      sgn = sgn && c != '-', c = getc();
    while (std::isdigit(c))
      x = x * 10 + c - '0', c = getc();
    x = sgn ? x : -x;
    return self;
  }
};

struct FastO : BasicBuffer {
  using BasicBuffer::BasicBuffer;
  std::array<char, 32> u{};
  ~FastO() {
    flush();
  }
  template <std::signed_integral T>
  FastO &operator<<(T x) {
    if (x < 0)
      putc('-'), x = -x;
    return self << std::make_unsigned(x);
  }
  template <std::unsigned_integral T>
  FastO &operator<<(T x) {
    char *i = u.data() + 30;
    do
      *--i = x % 10 + '0', x /= 10;
    while (x > 0);
    puts(i);
    return self;
  }
  FastO &operator<<(char x) {
    return putc(x), self;
  }
  FastO &operator<<(const char *x) {
    return puts(x), self;
  }
  FastO &operator<<(const std::string &x) {
    return puts(x.c_str()), self;
  }
};

#include "modint/modint-io.hpp"

#endif // ALGO_FASTIO
