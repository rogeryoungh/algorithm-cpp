#pragma once

#include "../base.hpp"
#include "../number/montgomery.hpp"
#include "../other/align-alloc.hpp"
#include "../number/barrett64.hpp"

#include <algorithm>
#include <span>

ALGO_BEGIN_NAMESPACE

template <class NTT>
struct CRT3Convolution {
  u32 P1 = 167772161, G1 = 3;
  u32 P2 = 469762049, G2 = 3;
  u32 P3 = 754974721, G3 = 11;
  Mont32 M1{P1}, M2{P2}, M3{P3};
  u32 ivm2 = M2.inv(M2.trans(P1));
  u32 ivm3 = M3.inv(M3.trans(P2));
  u32 ivm23 = M3.inv(M3.mul(M3.trans(P1), M3.trans(P2)));
  NTT ntt1{M1, G1}, ntt2{M2, G2}, ntt3{M3, G3};
  auto conv(std::span<const u32> f, std::span<const u32> g, const u32 mod) {

    u32 L = std::bit_ceil(f.size() + g.size() - 1);
    auto conv3 = [L, &f, &g](auto &&ntt) {
      AVec<u32> f1(L), g1(L);
      std::ranges::copy(f, f1.begin());
      std::ranges::copy(g, g1.begin());
      ntt.conv(f1.data(), g1.data(), L);
      return f1;
    };

    auto f1 = conv3(ntt1);
    auto f2 = conv3(ntt2);
    auto f3 = conv3(ntt3);

    Barrett64 bm(mod);
    const u64 M12 = u64(P1) * P2 % mod;
    for (u32 i = 0; i != L; ++i) {
      u32 a1 = f1[i], a2 = f2[i], a3 = f3[i];
      u32 x1 = M2.norm(M2.mul(M2.sub(a2, a1), ivm2));
      u32 x2 = M3.norm(M3.sub(M3.mul(ivm23, M3.sub(a3, a1)), M3.mul(x1, ivm3)));
      f1[i] = bm.calc(x2 * u64(M12) + x1 * u64(P1) + a1);
      // f1[i] = ((x2 * u64(P2) + x1) % mod * u64(P1) + a1) % mod;
    }
    f1.resize(f.size() + g.size() - 1);
    return f1;
  }
};

ALGO_END_NAMESPACE
