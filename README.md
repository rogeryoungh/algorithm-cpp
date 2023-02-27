# 算法库

[![verify](https://github.com/rogeryoungh/algorithm-cpp/actions/workflows/verify.yml/badge.svg)](https://github.com/rogeryoungh/algorithm-cpp/actions/workflows/verify.yml) [![code-size](https://img.shields.io/github/languages/code-size/rogeryoungh/algorithm-cpp)](https://github.com/rogeryoungh/algorithm-cpp)

> 正在编写中……这里是项目地址：[rogeryoungh/algorithm-cpp](https://github.com/rogeryoungh/algorithm-cpp)，感兴趣不妨来点个 star 哦！

我的算法库，计划包含算法竞赛的~~常见~~（我会的）算法。

本仓库的代码侧重于常数和易用性，封装较重。

如果你需要轻封装的代码，可以看看 [rogeryoungh/code-of-acm](https://github.com/rogeryoungh/code-of-acm)。

## 文档

还在施工。[文档](https://algo-cpp.rogery.dev/)。

## 使用方式

> TODO！

如果需要合并成但文件提交到 OJ，我提供了 `just` 的 `shell` 脚本，例如 `test/luogu/P3803.cpp`：

```bash
just expand test/luogu/P3803.cpp > P3803.full.cpp
```

已知问题：由于解析依赖 `gcc -E`，编译器会对所有 `include` 尝试展开，但又通过特殊方式排除了系统头文件，总之在 test 中引入系统头文件会导致错误。
