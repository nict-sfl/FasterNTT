//
// Created by Takuya HAYASHI on 2019/11/11.
//

#ifndef REDUCE_BENCH_H
#define REDUCE_BENCH_H

#include "cycle.h"
#include "gcc.h"
#include "x86_64.h"

#include <stdint.h>
#include <stdio.h>

#define OMITOPT(_out) { \
    FILE *_fp = fopen("/dev/null", "w");\
    fprintf(_fp, "%lu\n", _out);\
    fclose(_fp);\
}

void bench_mul128(const size_t niter) {
  const uint64_t p = p_kred;
  ticks t0, t1;
  uint64_t a, b, c;
  uint128_t t;
  const uint64_t *ptr = (uint64_t *) (&t);
  a = random_u64(p);
  b = random_u64(p);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
    t = mul128(a, b);
    c = ptr[0];
    t = mul128(b, c);
    a = ptr[0];
    t = mul128(c, a);
    b = ptr[0];
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_gcc(const size_t niter) {
  const uint64_t p = p_kred;
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p);
  b = random_u64(p);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
    c = mul_gcc(a, b, p);
    a = mul_gcc(b, c, p);
    b = mul_gcc(c, a, p);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_kred(const size_t niter) {
  const uint64_t p = p_kred, k = k_kred, s = s_kred;
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p);
  b = random_u64(p);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
    c = mul_kred(a, b, k, s);
    a = mul_kred(b, c, k, s);
    b = mul_kred(c, a, k, s);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_kred_asm(const size_t niter) {
  const uint64_t p = p_kred, k = k_kred, s = s_kred;
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p);
  b = random_u64(p);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
    c = mul_kred_asm(a, b, k, s);
    a = mul_kred_asm(b, c, k, s);
    b = mul_kred_asm(c, a, k, s);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_pm(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_pm);
  b = random_u64(p_pm);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
    c = mul_pm(a, b, k_pm, s_pm);
    a = mul_pm(b, c, k_pm, s_pm);
    b = mul_pm(c, a, k_pm, s_pm);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_pm2(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_pm);
  b = random_u64(p_pm);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
    c = mul_pm2(a, b, k_pm, s_pm);
    a = mul_pm2(b, c, k_pm, s_pm);
    b = mul_pm2(c, a, k_pm, s_pm);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_nfl(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_nfl);
  b = random_u64(p_nfl);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
    c = mul_nfl(a, b, p_nfl, v0, s0);
    a = mul_nfl(b, c, p_nfl, v0, s0);
    b = mul_nfl(c, a, p_nfl, v0, s0);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_nfl_asm(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_nfl);
  b = random_u64(p_nfl);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
    c = mul_nfl_asm(a, b, p_nfl, v0, s0);
    a = mul_nfl_asm(b, c, p_nfl, v0, s0);
    b = mul_nfl_asm(c, a, p_nfl, v0, s0);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_mm(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_mm);
  b = random_u64(p_mm);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
    c = mul_mm(a, b, p_mm, m_p_mm_i);
    a = mul_mm(b, c, p_mm, m_p_mm_i);
    b = mul_mm(c, a, p_mm, m_p_mm_i);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_mm_asm(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_mm);
  b = random_u64(p_mm);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
    c = mul_mm_asm(a, b, p_mm, m_p_mm_i);
    a = mul_mm_asm(b, c, p_mm, m_p_mm_i);
    b = mul_mm_asm(c, a, p_mm, m_p_mm_i);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_brt(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_brt);
  b = random_u64(p_brt);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
    c = mul_brt(a, b, p_brt, mu_brt);
    a = mul_brt(b, c, p_brt, mu_brt);
    b = mul_brt(c, a, p_brt, mu_brt);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_brt2(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c;
  a = random_u64(p_brt);
  b = random_u64(p_brt);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
    c = mul_brt2(a, b, p_brt, mu_brt2);
    a = mul_brt2(b, c, p_brt, mu_brt2);
    b = mul_brt2(c, a, p_brt, mu_brt2);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_brt_fixed(const size_t niter) {
  ticks t0, t1;
  uint64_t a, b, c, prep;
  a = random_u64(p_brt);
  b = random_u64(p_brt);
  prep = random_u64(p_brt);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
    c = mul_brt_fixed(a, b, p_brt, prep);
    a = mul_brt_fixed(b, c, p_brt, prep);
    b = mul_brt_fixed(c, a, p_brt, prep);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

void bench_mul_brt_fixed_asm(const size_t niter) {
  const uint64_t p = p_brt;
  ticks t0, t1;
  uint64_t a, b, c, bp;
  a = random_u64(p);
  b = random_u64(p);
  bp = random_u64(p);
  t0 = getticks();
  for (size_t i = 0; i < niter; ++i) {
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
    c = mul_brt_fixed_asm(a, b, bp, p);
    a = mul_brt_fixed_asm(b, c, bp, p);
    b = mul_brt_fixed_asm(c, a, bp, p);
  }
  t1 = getticks();
  OMITOPT(b);
  printf("%s: %.3lf\n", __FUNCTION__, elapsed(t1, t0) / 24 / niter);
}

#endif //REDUCE_BENCH_H
