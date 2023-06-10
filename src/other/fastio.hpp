#ifndef ALGO_FASTIO
#define ALGO_FASTIO

#include "../base.hpp"

#include <cstring>
#include <string>
#include <vector>

// ALGO_IO_NUMBER_ONLY 输入只有 0-9、-、空格、换行

namespace detail {

template <class Buf>
struct FastI : Buf {
  using Buf::pop;
  using Buf::top;
  FastI(FILE *f, u32 size = 1 << 18) : Buf(f, size) {}
  void skipSpace() {
#ifdef ALGO_IO_NUMBER_ONLY
    while (top() <= ' ')
      pop();
#else
    while (std::isspace(top()))
      pop();
#endif
  }
  FastI &operator>>(char &x) {
    skipSpace();
    x = pop();
    return *this;
  }
  FastI &operator>>(std::string &x) {
    x.resize(0);
    skipSpace();
    while (std::isgraph(top()))
      x.push_back(pop());
    return *this;
  }
  template <std::unsigned_integral T>
  FastI &operator>>(T &x) {
    x = 0;
    skipSpace();
#ifdef ALGO_IO_NUMBER_ONLY
    while (top() >= '0')
      x = x * 10 + (pop() & 0xf);
#else
    while (std::isdigit(top()))
      x = x * 10 + (pop() & 0xf);
#endif
    return *this;
  }
  template <std::signed_integral T>
  FastI &operator>>(T &x) {
    bool neg = false;
    x = 0;
    skipSpace();
    if (top() == '-')
      neg = true, pop();
#ifdef ALGO_IO_NUMBER_ONLY
    while (top() >= '0')
      x = x * 10 + (pop() & 0xf);
#else
    while (std::isdigit(top()))
      x = x * 10 + (pop() & 0xf);
#endif
    x = neg ? -x : x;
    return *this;
  }
};

template <class Buf>
struct FastO : Buf {
  using Buf::push;
  using Buf::push_uncheck;
  using Buf::puts;
  std::vector<u32> pre;
#define E10(x) u64(1E##x)
  FastO(FILE *f, u32 size = 1 << 18) : Buf(f, size), pre(E10(4)) {
    for (i32 i = 0; i < i32(E10(4)); ++i) {
      i32 ti = i;
      for (i32 j = 0; j < 4; ++j) {
        pre[i] = pre[i] << 8 | ti % 10 | 0x30;
        ti /= 10;
      }
    }
  }
  ~FastO() {
    Buf::flush();
  }
  template <std::signed_integral T>
  FastO &operator<<(T x) {
    if (x < 0)
      push('-'), x = -x;
    return *this << std::make_unsigned<T>::type(x);
  }
  void output4(u32 t) {
    auto tp = (const char *)&pre[t];
    if (t >= E10(2)) {
      if (t >= E10(3))
        push_uncheck(tp, 4);
      else
        push_uncheck(tp + 1, 3);
    } else {
      if (t >= E10(1))
        push_uncheck(tp + 2, 2);
      else
        push_uncheck(t | 0x30);
    }
  };
  template <std::unsigned_integral T>
  FastO &operator<<(T x) {
    Buf::reserve(32);
    if (x >= E10(8)) {
      u64 q0 = x / E10(8), r0 = x % E10(8);
      if (x >= E10(16)) {
        u64 q1 = q0 / E10(8), r1 = q0 % E10(8);
        output4(q1);
        push_uncheck(&pre[r1 / E10(4)], 4);
        push_uncheck(&pre[r1 % E10(4)], 4);
      } else if (x >= E10(12)) {
        output4(q0 / E10(4));
        push_uncheck(&pre[q0 % E10(4)], 4);
      } else {
        output4(q0);
      }
      push_uncheck(&pre[r0 / E10(4)], 4);
      push_uncheck(&pre[r0 % E10(4)], 4);
    } else {
      if (x >= E10(4)) {
        output4(x / E10(4));
        push_uncheck(&pre[x % E10(4)], 4);
      } else {
        output4(x);
      }
    }
    return *this;
  }
  FastO &operator<<(char x) {
    return push(x), *this;
  }
  FastO &operator<<(const char *x) {
    return puts(x), *this;
  }
  FastO &operator<<(const std::string &x) {
    return push(x.c_str(), x.size()), *this;
  }
};

struct BufO {
  FILE *f;
  char *beg, *end, *p;
  BufO(FILE *f_, u32 sz) : f(f_), beg(new char[sz]), end(beg + sz - 1), p(beg) {}
  ~BufO() {
    delete[] beg;
  }
  void flush() {
    std::fwrite(beg, 1, p - beg, f);
    p = beg;
  }
  void reserve(u32 len) {
    if (end - p <= i32(len))
      flush();
  }
  void push(char s) {
    *p++ = s;
    reserve(0);
  }
  void push(const char *s, u32 len) {
    reserve(len);
    push_uncheck(s, len);
  }
  void push_uncheck(char s) {
    *p++ = s;
  }
  void push_uncheck(const void *s, u32 len) {
    std::memcpy(p, s, len);
    p += len;
  }
  void puts(const char *s) {
    while (*s != 0)
      push(*s++);
  }
};

} // namespace detail

#if defined(_WIN32) || defined(WIN32)

namespace detail {

struct BufI {
  FILE *f;
  char *beg, *end, *p;
  BufI(FILE *f_, u32 sz) : f(f_), beg(new char[sz]), end(beg + sz - 1), p(end) {
    reread();
  }
  ~BufI() {
    delete[] beg;
  }
  void reread() {
    std::memmove(beg, p, end - p);
    p = end - p + beg;
    u32 r = std::fread(p, 1, end - p, f);
    p[r] = 0;
    p = beg;
  }
  char pop() {
    char r = *p++;
    if (p == end)
      reread();
    return r;
  }
  char top() const {
    return *p;
  }
};

} // namespace detail

#else

#include <sys/mman.h>
#include <sys/stat.h>

namespace detail {

struct BufI {
  struct stat sb;
  char *p;
  BufI(FILE *f, u32) {
    i32 fd = fileno(f);
    fstat(fd, &sb);
    p = (char *)mmap(nullptr, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    madvise(p, sb.st_size, MADV_SEQUENTIAL);
  }
  ~BufI() {
    munmap(p, sb.st_size);
  }
  char pop() {
    return *p++;
  }
  char top() const {
    return *p;
  }
};

} // namespace detail

#endif

using FastI = detail::FastI<detail::BufI>;

using FastO = detail::FastO<detail::BufO>;

#include "modint/modint-io.hpp"

#endif // ALGO_FASTIO
