#ifndef ALGO_H_FASTIO
#define ALGO_H_FASTIO

#include "../base.hpp"

#include <array>
#include <bit>
#include <cstring>
#include <string>
#include <vector>

// ALGO_IO_NUMBER_ONLY 输入只有 0-9、-、空格、换行

#if defined(_WIN32) || defined(WIN32)

ALGO_BEGIN_NAMESPACE

struct FastI {
  FILE *const f;
  const usize bufsz;
  std::vector<char> buf;
  usize p{};
  FastI(FILE *file, usize size = 1 << 18) : f(file), bufsz(size), buf(size + 4), p(size) {
    reread();
  }
  void reserve(usize len) {
    if (bufsz - p < len)
      reread();
  }
  char pop() {
    reserve(4);
    return buf[p++];
  }
  char top() const {
    return buf[p];
  }
  u32 top4() {
    u32 t;
    std::memcpy(&t, &buf[p], sizeof(t));
    return t;
  }
  void reread() {
    std::memmove(&buf[0], &buf[p], bufsz - p);
    u32 cnt = std::fread(&buf[bufsz - p], 1, p, f);
    buf[bufsz - p + cnt] = '\0', p = 0;
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
    p = (char *)mmap(nullptr, sb.st_size + 4, PROT_READ, MAP_PRIVATE, fd, 0);
    madvise(p, sb.st_size, MADV_SEQUENTIAL);
  }
  ~FastI() {
    munmap(p, sb.st_size + 4);
  }
  char pop() {
    return *p++;
  }
  char top() {
    return *p;
  }
  void reserve(usize) const {}
  u32 top4() {
    u32 t;
    std::memcpy(&t, p, sizeof(t));
    return t;
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
  FastO(FILE *file, u32 size = 1 << 18) : f(file), bufsz(size), buf(size + 4) {
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
      reserve(len);
      std::memcpy(&buf[p], s, len);
      p += len;
    } else {
      flush();
      std::fwrite(s, 1, len, f);
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
    constexpr std::array<u32, 5> p10 = {u64(1E0), u64(1E1), u64(1E2), u64(1E3), u64(1E4)};
    while (true) {
      fin.reserve(4);
      u32 u = fin.top4();
#ifdef ALGO_IO_NUMBER_ONLY
      u32 umask = u & 0xf0f0f0f0;
#else
      u32 umask = u & (u + 0x06060606) & 0xf0f0f0f0;
#endif
      u32 len = std::countr_zero(umask ^ 0x30303030) >> 3;
      if (len == 0)
        break;
      u <<= sizeof(u) * 8 - (len << 3);
      u = (u & 0x0f0f0f0f) * 0x0a01 >> 0x08;
      u = (u & 0x00ff00ff) * 0x00640001 >> 0x10;
      fin.p += len;
      x = x * p10[len] + u;
      if (len != sizeof(u))
        break;
    }
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
