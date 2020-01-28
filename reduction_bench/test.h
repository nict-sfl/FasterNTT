//
// Created by Takuya HAYASHI on 2019/11/11.
//

#ifndef REDUCE_TEST_H
#define REDUCE_TEST_H

#include "gcc.h"
#include "x86_64.h"

#include <stdint.h>
#include <stdio.h>

void test_pm(const size_t niter) {
  const uint64_t p = p_pm, k = k_pm, s = s_pm;
  uint64_t a, b, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_pm(a, b, k, s);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_pm2(const size_t niter) {
  const uint64_t p = p_pm, k = k_pm, s = s_pm;
  uint64_t a, b, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_pm2(a, b, k, s);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}


void test_mm(const size_t niter) {
  const uint64_t p = p_mm, m = m_p_mm_i;
  uint64_t a, a_, b, b_, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    a_ = mul_mm(a, rr, p, m);
    b_ = mul_mm(b, rr, p, m);
    c_ = mul_mm(a_, b_, p, m);
    c_ = mul_mm(c_, 1, p, m);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}


void test_mm_asm(const size_t niter) {
  const uint64_t p = p_mm, m = m_p_mm_i;
  uint64_t a, a_, b, b_, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    a_ = mul_mm_asm(a, rr, p, m);
    b_ = mul_mm_asm(b, rr, p, m);
    c_ = mul_mm_asm(a_, b_, p, m);
    c_ = mul_mm_asm(c_, 1, p, m);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_brt(const size_t niter) {
  const uint64_t p = p_brt, *mu = mu_brt;
  uint64_t a, b, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_brt(a, b, p, mu);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_brt2(const size_t niter) {
  const uint64_t p = p_brt2, mu = mu_brt2;
  uint64_t a, b, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_brt2(a, b, p, mu);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_brt_fixed(const size_t niter) {
  const uint64_t p = p_brt;
  uint64_t a, b, bp, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);
    bp = mul_brt_prep(b, p);

    c = mul_gcc(a, b, p);

    c_ = mul_brt_fixed(a, b, bp, p);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_brt_fixed_asm(const size_t niter) {
  const uint64_t p = p_brt;
  uint64_t a, b, bp, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);
    bp = mul_brt_prep(b, p);

    c = mul_gcc(a, b, p);

    c_ = mul_brt_fixed_asm(a, b, bp, p);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_nfl(const size_t niter) {
  const uint64_t p = p_nfl, s = s0, v = v0;
  uint64_t a, b, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_nfl(a, b, p, v0, s0);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_nfl_asm(const size_t niter) {
  static const uint64_t p = p_nfl, s = s0, v = v0;
  uint64_t a, b, c, c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_nfl_asm(a, b, p, v0, s0);
#ifndef REDC_LAZY
    if(c != c_){
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#else
    if (c != c_ % p) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
#endif
  }
  printf("%s: o\n", __FUNCTION__);
}

uint64_t corr_kred(const int64_t c, const uint64_t k_, const uint64_t p) {
  uint128_t t = mul128((c % (const int64_t) p) + p, k_);
  t.d128 %= p;
  return t.d64[0];
}

void test_kred(const size_t niter) {
  const uint64_t p = p_kred, s = s_kred, k = k_kred, k_ = k__kred;
  uint64_t a, b, c;
  int64_t c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_kred(a, b, k, s);
    c_ = corr_kred(c_, k_, p);
    if (c != c_) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
  }
  printf("%s: o\n", __FUNCTION__);
}

void test_kred_asm(const size_t niter) {
  const uint64_t p = p_kred, s = s_kred, k = k_kred, k_ = k__kred;
  uint64_t a, b, c;
  int64_t c_;
  for (size_t i = 0; i < niter; ++i) {
    a = random_u64(p);
    b = random_u64(p);

    c = mul_gcc(a, b, p);

    c_ = mul_kred_asm(a, b, k, s);
    c_ = corr_kred(c_, k_, p);
    if (c != c_) {
      printf("%s: x... %lu: %lu %lu %lu %lu\n", __FUNCTION__, i, a, b, c, c_);
      return;
    }
  }
  printf("%s: o\n", __FUNCTION__);
}

#endif //REDUCE_TEST_H
