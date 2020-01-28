//
// Created by Takuya HAYASHI on 2019/11/07.
//

#include <memory>
#include <benchmark/benchmark.h>
#include <yell/poly.hpp>
#include <nfl/poly.hpp>
#include "Rq.hpp"

static constexpr size_t qNum = 4;
static constexpr size_t N = 8192;
static constexpr size_t NRep = 64;

fntt::Rq<N, qNum> a1, b1, c1;
nfl::poly<uint64_t, N, qNum> a2, b2, c2;
yell::poly<N> a3 = yell::poly<N>(qNum, yell::uniform()),
  b3 = yell::poly<N>(qNum, yell::uniform()),
  c3 = yell::poly<N>(qNum, yell::uniform());

struct NTT : public benchmark::Fixture {
  void SetUp(const ::benchmark::State &st) {
    yell::ntt<N>::init_ntt_tables(qNum);
    a1.rand();
    b1.rand();
    for(int cm = 0; cm < qNum; ++cm) {
      std::memcpy(a2.begin() + N * cm, a3.cptr_at(cm), sizeof(uint64_t) * N);
      std::memcpy(b2.begin() + N * cm, b3.cptr_at(cm), sizeof(uint64_t) * N);
    }
  }

  void TearDown(const ::benchmark::State &st) {
  }
};

BENCHMARK_DEFINE_F(NTT, FasterNTT_NTT)(benchmark::State &st) {
  for(auto _ : st) a1.ntt();
}

BENCHMARK_DEFINE_F(NTT, FasterNTT_INTT)(benchmark::State &st) {
  for(auto _ : st) b1.intt();
}

BENCHMARK_DEFINE_F(NTT, FasterNTT_PMUL)(benchmark::State &st) {
  for(auto _ : st) fntt::Rq<N, qNum>::pmul_lazy(c1, a1, b1);
}

BENCHMARK_DEFINE_F(NTT, NFLlib_NTT)(benchmark::State &st) {
  for(auto _ : st) a2.ntt_pow_phi();
}

BENCHMARK_DEFINE_F(NTT, NFLlib_INTT)(benchmark::State &st) {
  for(auto _ : st) b2.invntt_pow_invphi();
}

BENCHMARK_DEFINE_F(NTT, NFLlib_PMUL)(benchmark::State &st) {
  for(auto _ : st) c2 = a2 * b2;
}

BENCHMARK_DEFINE_F(NTT, YELL_NTT)(benchmark::State &st) {
  for(auto _ : st) a3.forward_lazy();
}

BENCHMARK_DEFINE_F(NTT, YELL_INTT)(benchmark::State &st) {
  for(auto _ : st) b3.backward();
}

BENCHMARK_DEFINE_F(NTT, YELL_PMUL)(benchmark::State &st) {
  for(auto _ : st) c3 = a3 * b3;
}

BENCHMARK_REGISTER_F(NTT, FasterNTT_NTT)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FasterNTT_INTT)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FasterNTT_PMUL)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, NFLlib_NTT)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, NFLlib_INTT)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, NFLlib_PMUL)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, YELL_NTT)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, YELL_INTT)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, YELL_PMUL)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_MAIN();
