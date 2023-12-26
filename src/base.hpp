#ifndef ALGO_H_BASE
#define ALGO_H_BASE

// GENERATE FROM: https://github.com/rogeryoungh/algorithm-cpp

#ifndef ALGO_NO_NAMESPACE
#define ALGO_BEGIN_NAMESPACE namespace algo {
#define ALGO_END_NAMESPACE }
#else
#define ALGO_BEGIN_NAMESPACE
#define ALGO_END_NAMESPACE
#endif

#include <cstdint>

ALGO_BEGIN_NAMESPACE

using i8 = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;
using i128 = __int128_t;

using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

using usize = std::size_t;

using f32 = float;
using f64 = double;
using f80 = long double;

ALGO_END_NAMESPACE

#endif
