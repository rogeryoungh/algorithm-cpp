#ifndef ALGO_ALIGN_ALLOC
#define ALGO_ALIGN_ALLOC

#include "../base.hpp"

#include <mm_malloc.h>
#include <vector>

ALGO_BEGIN_NAMESPACE

template <class T, u32 align>
struct AlignedAllocator {
  using value_type = T;
  static T *allocate(std::size_t n) {
    // return new (std::align_val_t(align)) T[n];
    return reinterpret_cast<T *>(std::aligned_alloc(align, n * sizeof(T)));
  }
  template <class U>
  struct rebind {
    using other = AlignedAllocator<U, align>;
  };
  static void deallocate(T *p, std::size_t) {
    // ::operator delete[](p, std::align_val_t(align));
    std::free(p);
  }
};

template <class T, u32 align = 32>
using AVec = std::vector<T, AlignedAllocator<T, align>>;

ALGO_END_NAMESPACE

#endif // ALGO_ALIGN_ALLOC
