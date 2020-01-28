/*
 * @file test.cpp 
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 */

#include "Rq.hpp"

using namespace fntt;

#include <NTL/ZZ_pX.h>
template<size_t N, size_t qNum, enum NTTSelector NTT, enum ImplSelector Impl>
void MulNTL(Rq<N, qNum, NTT, Impl> &c, const Rq<N, qNum, NTT, Impl> &a, const Rq<N, qNum, NTT, Impl> &b) {
  for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
    NTL::ZZ_p::init(NTL::conv<NTL::ZZ>((unsigned long)Rq<N, qNum, NTT, Impl>::param_.p[qIdx]));
    NTL::ZZ_pX f;
    f.SetLength(N + 1);
    f[N] = f[0] = 1;
    NTL::ZZ_pXModulus F(f);

    NTL::ZZ_pX a_, b_, c_;
    a_.SetLength(N);
    b_.SetLength(N);
    c_.SetLength(N);

    for(size_t i = 0; i < N; ++i) a_[i] = a(qIdx, i);
    for(size_t i = 0; i < N; ++i) b_[i] = b(qIdx, i);

    NTL::MulMod(c_, a_, b_, F);

    for(size_t i = 0; i < N; ++i) c(qIdx, i) = NTL::conv<long>(c_[i]);
  }
}

template <size_t N, size_t qNum, enum NTTSelector NTT, enum ImplSelector Impl>
bool testMulMax() {
  using R = Rq<N, qNum, NTT, Impl>;
  R a, b, c, c_;
  for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
    for(size_t i = 0; i < N; ++i) {
      a(qIdx, i) = R::param_.p[qIdx] - 1UL;
      b(qIdx, i) = R::param_.p[qIdx] - 1UL;
    }
  }

  MulNTL<N, qNum, NTT, Impl>(c, a, b);
  R::mul(c_, a, b);

  return c != c_;
}

template<size_t N, size_t qNum, enum NTTSelector NTT, enum ImplSelector Impl>
void TestMul(const size_t niter = 128) {
  using R = Rq<N, qNum, NTT, Impl>;
  R a, b, c, c_;

  printf("mul(%lu %lu %d %d): ", N, qNum, NTT, Impl);

  if(testMulMax<N, qNum, NTT, Impl>()) {
    printf("x max\n");
    return;
  }

  for(size_t i = 0; i < niter; ++i) {
    a.rand();
    b.rand();

    MulNTL<N, qNum, NTT, Impl>(c, a, b);
    R::mul(c_, a, b);

    if(c != c_) {
      printf("x rand\n");
      return;
    }
  }
  printf("o\n");
}

static constexpr size_t qNum = 2;

int main() {

  TestMul<4096, qNum, LN16, Generic>();
  TestMul<8192, qNum, LN16, Generic>();
  TestMul<16384, qNum, LN16, Generic>();
  TestMul<32768, qNum, LN16, Generic>();

  TestMul<4096, qNum, LN16, x86_64>();
  TestMul<8192, qNum, LN16, x86_64>();
  TestMul<16384, qNum, LN16, x86_64>();
  TestMul<32768, qNum, LN16, x86_64>();

  TestMul<4096, qNum, S17, Generic>();
  TestMul<8192, qNum, S17, Generic>();
  TestMul<16384, qNum, S17, Generic>();
  TestMul<32768, qNum, S17, Generic>();

  TestMul<4096, qNum, S17, x86_64>();
  TestMul<8192, qNum, S17, x86_64>();
  TestMul<16384, qNum, S17, x86_64>();
  TestMul<32768, qNum, S17, x86_64>();

  TestMul<4096, qNum, FNTT, Generic>();
  TestMul<8192, qNum, FNTT, Generic>();
  TestMul<16384, qNum, FNTT, Generic>();
  TestMul<32768, qNum, FNTT, Generic>();

  TestMul<4096, qNum, FNTT, x86_64>();
  TestMul<8192, qNum, FNTT, x86_64>();
  TestMul<16384, qNum, FNTT, x86_64>();
  TestMul<32768, qNum, FNTT, x86_64>();

  TestMul<4096, qNum, FNTT, SSE>();
  TestMul<8192, qNum, FNTT, SSE>();
  TestMul<16384, qNum, FNTT, SSE>();
  TestMul<32768, qNum, FNTT, SSE>();

  TestMul<4096, qNum, FNTT, AVX2>();
  TestMul<8192, qNum, FNTT, AVX2>();
  TestMul<16384, qNum, FNTT, AVX2>();
  TestMul<32768, qNum, FNTT, AVX2>();

  return 0;
}