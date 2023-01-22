---
author: "rogeryoungh"
title: "Baillie PSW 素性测试"
date: "2023-01-22"
---

## 理论

理论我也讲不出啥。

一个显然的观察是，Fermat 反素数和 Lucas 反素数的从计算过程来看就很不一样，换基不如换算法。把两个算法糅到一块就行了。



## 算法描述

输入整数 $N$。

1. 试除小素数，若能整除则输出 $N$ 为合数。
2. 使用基为 $2$ 的 Miller Rabin 算法
3. 进行一次 Lucas 素性测试判断。

$2^{64}$ 内无反例。

## 备注

在 LOJ143 上跑起来没 Miller Rabin 算法快，也可能是数据有所针对。

## 参考

- [素数判定 - mathmu](http://mathmu.github.io/MTCAS/doc/PrimeTest.html)
- [Baillie–PSW primality test - Wikipedia](https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test)

## 代码

{{<code file="./baillie-psw.hpp" lang="cpp">}}
