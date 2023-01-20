/////////////////////////////////////////////////////////

#include <benchmark/benchmark.h>
#include <chrono>

struct TEST_BASE {
  void init(int n);
  int run(int n);
};

template <class T>
class BM_TEST_F : public benchmark::Fixture {
public:
  void SetUp(const benchmark::State &state) override {
    runner.init(state.range());
  }
  void test(benchmark::State &state) {
    state.counters["ntt"] =
        benchmark::Counter(runner.run(state.range()), benchmark::Counter::kDefaults,
                           benchmark::Counter::OneK::kIs1024);
  }
  T runner;
};

#define BM_DEF(X)                                                                 \
  BENCHMARK_TEMPLATE_DEFINE_F(BM_TEST_F, X, X)(benchmark::State & state) {        \
    for (auto _ : state) {                                                        \
      auto start = std::chrono::high_resolution_clock::now();                     \
      test(state);                                                                \
      auto end = std::chrono::high_resolution_clock::now();                       \
      auto elapsed_seconds =                                                      \
          std::chrono::duration_cast<std::chrono::duration<double>>(end - start); \
    }                                                                             \
                                                                                  \
    state.SetComplexityN(state.range());                                          \
  }                                                                               \
  BENCHMARK_REGISTER_F(BM_TEST_F, X)->Unit(benchmark::kMillisecond)->Complexity()->Name(#X)

BENCHMARK_MAIN();
