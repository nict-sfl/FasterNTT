/*
 * @file bench.cpp 
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 */

#include<benchmark/benchmark.h>
#include "Rq.hpp"

using namespace fntt;

struct FNTT : public benchmark::Fixture {
  static constexpr unsigned int N = 32768;
  using R = Rq<N, 4, fntt::FNTT, fntt::AVX2>;
  R a, b, c;

  void SetUp(const ::benchmark::State &state) {
    a.rand();
    b.rand();
  }

  void TearDown(const ::benchmark::State &state) {
  }
};

BENCHMARK_DEFINE_F(FNTT, ntt)(benchmark::State &st) {
  for (auto _ : st) a.ntt();
}

BENCHMARK_DEFINE_F(FNTT, intt)(benchmark::State &st) {
  for (auto _ : st) a.intt();
}

BENCHMARK_DEFINE_F(FNTT, pmul)(benchmark::State &st) {
  for (auto _ : st) R::pmul(c, a, b);
}

constexpr size_t NRep = 64;

BENCHMARK_REGISTER_F(FNTT, ntt)->Unit(benchmark::kMicrosecond)->Threads(1)->
    Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(FNTT, intt)->Unit(benchmark::kMicrosecond)->Threads(1)->
    Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(FNTT, pmul)->Unit(benchmark::kMicrosecond)->Threads(1)->
    Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_MAIN();
