#ifndef ALGO_H_FASTIO
#define ALGO_H_FASTIO

#include "../base.hpp"

#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <bit>

// ALGO_IO_NUMBER_ONLY 输入只有 0-9、-、空格、换行

#if defined(_WIN32) || defined(WIN32)

ALGO_BEGIN_NAMESPACE

struct FastI {
  FILE *const f;
  const usize bufsz;
  std::vector<char> buf;
  usize p{};
  FastI(FILE *file, usize size = 1 << 18) : f(file), bufsz(size - 1), buf(size), p(bufsz) {
    reread();
  }
  char pop() {
    char r = buf[p++];
    if (p == bufsz)
      reread();
    return r;
  }
  void reread() {
    std::memmove(&buf[0], &buf[p], bufsz - p);
    u32 cnt = std::fread(&buf[bufsz - p], 1, p, f);
    buf[cnt] = '\0', p = 0;
  }
  bool static isspace(char c) {
#ifdef ALGO_IO_NUMBER_ONLY
    return c <= ' ';
#else
    return std::isspace(c);
#endif
  }
  i32 skipSpace() {
    i32 r = pop();
    while (isspace(r))
      r = pop();
    return r;
  }
};

#else

#include <sys/mman.h>
#include <sys/stat.h>

ALGO_BEGIN_NAMESPACE

struct FastI {
  struct stat sb;
  char *p;
  FastI(FILE *file) {
    i32 fd = fileno(file);
    fstat(fd, &sb);
    p = (char *)mmap(nullptr, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    madvise(p, sb.st_size, MADV_SEQUENTIAL);
  }
  ~FastI() {
    munmap(p, sb.st_size);
  }
  char pop() {
    return *p++;
  }
  char top() {
    return *p;
  }
  void skipSpace() {
#ifdef ALGO_IO_NUMBER_ONLY
    while (top() <= ' ')
      pop();
#else
    while (std::isspace(top()))
      pop();
#endif
  }
};

#endif

struct FastO {
  FILE *const f;
  const usize bufsz;
  std::vector<char> buf;
  std::array<u32, 10000> pre;
  usize p{};
  FastO(FILE *file, u32 size = 1 << 18) : f(file), bufsz(size - 1), buf(size) {
    for (u32 i = 0; i != 10000; ++i) {
      u32 k = i;
      for (u32 j = 0; j != 4; ++j) {
        pre[i] = pre[i] << 8 | k % 10 | '0';
        k /= 10;
      }
    }
  }
  ~FastO() {
    flush();
  }
  void flush() {
    std::fwrite(&buf[0], 1, p, f);
    p = 0;
  }
  void reserve(usize len) {
    if (bufsz - p < len)
      flush();
  }
  void push(const char *s, usize len) {
    if (len < bufsz / 2) {
      reserve(len + 1);
      std::memcpy(&buf[p], s, len);
      p += len;
    } else {
      flush();
      fwrite(s, 1, len, f);
    }
  }
  template <u32 check = 1>
  void push(char c) {
    if constexpr (check)
      reserve(1);
    buf[p++] = c;
  }
};

template <class T>
inline FastI &operator>>(FastI &fin, T &x) {
  if constexpr (std::is_same_v<T, char>) {
    fin.skipSpace();
    x = fin.pop();
  } else if constexpr (std::is_same_v<T, std::string>) {
    x.resize(0);
    while (std::isgraph(fin.top()))
      x.push_back(fin.pop());
  } else if constexpr (std::is_integral_v<T>) {
    x = 0;
    bool neg = false;
    fin.skipSpace();
    if (std::is_signed_v<T> && fin.top() == '-')
      fin.pop(), neg = true;
#ifdef ALGO_IO_NUMBER_ONLY
    constexpr std::array<u64, 9> p10 = {u64(1E0), u64(1E1), u64(1E2), u64(1E3), u64(1E4),
                                        u64(1E5), u64(1E6), u64(1E7), u64(1E8)};
    while (true) {
      u64 u = 0;
      std::memcpy(&u, fin.p, sizeof(u));
      u64 len = std::countr_zero((u & 0xf0f0f0f0f0f0f0f0) ^ 0x3030303030303030) >> 3;
      if (len == 0)
        break;
      u <<= sizeof(u) * 8 - (len << 3);
      u = (u & 0x0f0f0f0f0f0f0f0f) * 0x0a01 >> 0x08;
      u = (u & 0x00ff00ff00ff00ff) * 0x00640001 >> 0x10;
      u = (u & 0x0000ffff0000ffff) * 0x271000000001 >> 0x20;
      fin.p += len;
      x = x * p10[len] + u;
      if (len != sizeof(u))
        break;
    }
#else
    while (std::isdigit(fin.top()))
      x = x * 10 + (fin.pop() & 0xf);
#endif
    if constexpr (std::is_signed_v<T>)
      x = neg ? -x : x;
  }
  return fin;
}

template <class T>
inline FastO &operator<<(FastO &fout, const T &x) {
  if constexpr (std::is_array_v<T>) {
    fout.push(x, sizeof(x) - 1);
  } else if constexpr (std::is_same_v<T, char>) {
    fout.push(x);
  } else if constexpr (std::is_integral_v<T>) {
    if (x == 0) {
      fout.push('0');
    } else {
      std::make_unsigned_t<T> nx = x;
      if (std::is_signed_v<T> && x < 0)
        fout.push('-'), nx = -nx;
      std::array<char, 32> arr{};
      usize s = 32;
      while (nx >= 10000) {
        std::memcpy(&arr[s -= 4], &fout.pre[nx % 10000], 4);
        nx /= 10000;
      }
      std::memcpy(&arr[s -= 4], &fout.pre[nx], 4);
      s += nx < 100 ? (nx < 10 ? 3 : 2) : (nx < 1000 ? 1 : 0);
      fout.push(&arr[s], 32 - s);
    }
  }
  return fout;
}

ALGO_END_NAMESPACE

#endif
