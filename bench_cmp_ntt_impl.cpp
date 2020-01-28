/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#include <memory>
#include <benchmark/benchmark.h>
#include "Rq.hpp"

static constexpr size_t qNum = 4;
static constexpr size_t N = 8192;
static constexpr size_t NRep = 64;

using R1 = fntt::Rq<N, qNum, fntt::LN16, fntt::Generic>;
using R2 = fntt::Rq<N, qNum, fntt::S17, fntt::Generic>;
using R3 = fntt::Rq<N, qNum, fntt::FNTT, fntt::Generic>;
using R4 = fntt::Rq<N, qNum, fntt::LN16, fntt::x86_64>;
using R5 = fntt::Rq<N, qNum, fntt::S17, fntt::x86_64>;
using R6 = fntt::Rq<N, qNum, fntt::FNTT, fntt::x86_64>;
using R7 = fntt::Rq<N, qNum, fntt::FNTT, fntt::SSE>;
using R8 = fntt::Rq<N, qNum, fntt::FNTT, fntt::AVX2>;

R1 a1, b1, c1;
R2 a2, b2, c2;
R3 a3, b3, c3;
R4 a4, b4, c4;
R5 a5, b5, c5;
R6 a6, b6, c6;
R7 a7, b7, c7;
R8 a8, b8, c8;

struct NTT : public benchmark::Fixture {
  void SetUp(const ::benchmark::State &st) {
    a1.rand(), b1.rand();
    a2.rand(), b2.rand();
    a3.rand(), b3.rand();
    a4.rand(), b4.rand();
    a5.rand(), b5.rand();
    a6.rand(), b6.rand();
    a7.rand(), b7.rand();
    a8.rand(), b8.rand();
  }

  void TearDown(const ::benchmark::State &st) {
  }
};

BENCHMARK_DEFINE_F(NTT, LongaNaehrig16_NTT_Generic)(benchmark::State &st) {
  for(auto _ : st) a1.ntt();
}

BENCHMARK_DEFINE_F(NTT, LongaNaehrig16_INTT_Generic)(benchmark::State &st) {
  for(auto _ : st) a1.intt();
}

BENCHMARK_DEFINE_F(NTT, LongaNaehrig16_PMUL_Generic)(benchmark::State &st) {
  for(auto _ : st) R1::pmul_lazy(c1, a1, b1);
}

BENCHMARK_DEFINE_F(NTT, Scott17_NTT_Generic)(benchmark::State &st) {
  for(auto _ : st) a2.ntt();
}

BENCHMARK_DEFINE_F(NTT, Scott17_INTT_Generic)(benchmark::State &st) {
  for(auto _ : st) a2.intt();
}

BENCHMARK_DEFINE_F(NTT, Scott17_PMUL_Generic)(benchmark::State &st) {
  for(auto _ : st) R2::pmul_lazy(c2, a2, b2);
}

BENCHMARK_DEFINE_F(NTT, FNTT_NTT_Generic)(benchmark::State &st) {
  for(auto _ : st) a3.ntt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_INTT_Generic)(benchmark::State &st) {
  for(auto _ : st) a3.intt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_PMUL_Generic)(benchmark::State &st) {
  for(auto _ : st) R3::pmul_lazy(c3, a3, b3);
}

BENCHMARK_DEFINE_F(NTT, LongaNaehrig16_NTT_x86_64)(benchmark::State &st) {
  for(auto _ : st) a4.ntt();
}

BENCHMARK_DEFINE_F(NTT, LongaNaehrig16_INTT_x86_64)(benchmark::State &st) {
  for(auto _ : st) a4.intt();
}

BENCHMARK_DEFINE_F(NTT, LongaNaehrig16_PMUL_x86_64)(benchmark::State &st) {
  for(auto _ : st) R4::pmul_lazy(c4, a4, b4);
}

BENCHMARK_DEFINE_F(NTT, Scott17_NTT_x86_64)(benchmark::State &st) {
  for(auto _ : st) a5.ntt();
}

BENCHMARK_DEFINE_F(NTT, Scott17_INTT_x86_64)(benchmark::State &st) {
  for(auto _ : st) a5.intt();
}

BENCHMARK_DEFINE_F(NTT, Scott17_PMUL_x86_64)(benchmark::State &st) {
  for(auto _ : st) R5::pmul_lazy(c5, a5, b5);
}

BENCHMARK_DEFINE_F(NTT, FNTT_NTT_x86_64)(benchmark::State &st) {
  for(auto _ : st) a6.ntt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_INTT_x86_64)(benchmark::State &st) {
  for(auto _ : st) a6.intt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_PMUL_x86_64)(benchmark::State &st) {
  for(auto _ : st) R6::pmul_lazy(c6, a6, b6);
}

BENCHMARK_DEFINE_F(NTT, FNTT_NTT_SSE)(benchmark::State &st) {
  for(auto _ : st) a7.ntt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_INTT_SSE)(benchmark::State &st) {
  for(auto _ : st) a7.intt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_PMUL_SSE)(benchmark::State &st) {
  for(auto _ : st) R7::pmul_lazy(c7, a7, b7);
}

BENCHMARK_DEFINE_F(NTT, FNTT_NTT_AVX2)(benchmark::State &st) {
  for(auto _ : st) a8.ntt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_INTT_AVX2)(benchmark::State &st) {
  for(auto _ : st) a8.intt();
}

BENCHMARK_DEFINE_F(NTT, FNTT_PMUL_AVX2)(benchmark::State &st) {
  for(auto _ : st) R8::pmul_lazy(c8, a8, b8);
}

BENCHMARK_REGISTER_F(NTT, LongaNaehrig16_NTT_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, LongaNaehrig16_INTT_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, LongaNaehrig16_PMUL_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, Scott17_NTT_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, Scott17_INTT_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, Scott17_PMUL_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, FNTT_NTT_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_INTT_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_PMUL_Generic)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, LongaNaehrig16_NTT_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, LongaNaehrig16_INTT_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, LongaNaehrig16_PMUL_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, Scott17_NTT_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, Scott17_INTT_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, Scott17_PMUL_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, FNTT_NTT_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_INTT_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_PMUL_x86_64)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, FNTT_NTT_SSE)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_INTT_SSE)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_PMUL_SSE)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_REGISTER_F(NTT, FNTT_NTT_AVX2)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_INTT_AVX2)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);
BENCHMARK_REGISTER_F(NTT, FNTT_PMUL_AVX2)->Unit(benchmark::kMicrosecond)->Threads(1)->
  Repetitions(NRep)->DisplayAggregatesOnly(true);

BENCHMARK_MAIN();
