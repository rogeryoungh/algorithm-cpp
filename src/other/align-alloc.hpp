#ifndef ALGO_ALIGN_ALLOC
#define ALGO_ALIGN_ALLOC

#include "../base.hpp"

#include <cstdlib>
#include <vector>
#include <new>

template <class T, u32 align>
struct AlignedAllocator {
  using value_type = T;
  T *allocate(std::size_t n) {
    return new (std::align_val_t(align)) T[n];
  }
  template <class U>
  struct rebind {
    using other = AlignedAllocator<U, align>;
  };
  void deallocate(T *p, std::size_t) {
    ::operator delete[] (p, std::align_val_t(align));
  }
};

template <class T, u32 align = 32>
using AVec = std::vector<T, AlignedAllocator<T, align>>;

#endif // ALGO_ALIGN_ALLOC
