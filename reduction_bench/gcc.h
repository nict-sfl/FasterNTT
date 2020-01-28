//
// Created by Takuya HAYASHI on 2019/11/11.
//

#ifndef REDUCE_GCC_H
#define REDUCE_GCC_H

#include "param.h"
#include <gmp.h>

typedef union uint128_data {
  unsigned __int128 d128;
  uint64_t d64[2];
} uint128_t;

uint64_t random_u64(const uint64_t p){
  uint64_t r;
  mpn_random((long unsigned int*)&r, 1);
  return r % p;
}

inline uint128_t mul128(const uint64_t a, const uint64_t b){
  uint128_t r;
  r.d128 = (unsigned __int128)a * b;
  return r;
}

inline uint128_t lshift128(uint128_t a, const size_t s){
  a.d128 <<= s;
  return a;
}

inline uint128_t rshift128(uint128_t a, const size_t s){
  a.d128 >>= s;
  return a;
}

inline uint128_t adds128_64(uint128_t a, const uint64_t b){
  a.d128 += b;
  return a;
}

inline uint128_t adds128(uint128_t a, const uint128_t b){
  a.d128 += b.d128;
  return a;
}

inline uint128_t subs128(uint128_t a, const uint128_t b){
  a.d128 -= b.d128;
  return a;
}

uint64_t mul_gcc(const uint64_t a, const uint64_t b, const uint64_t p){
  return ((unsigned __int128)a * b) % p;
}

int64_t kred(const uint128_t c, const uint64_t k, const size_t s) {
  uint64_t t0, t1;
  t0 = c.d64[0] & ((1UL << s) - 1UL);
  t1 = (c.d64[1] << (WORD_SIZE - s)) | (c.d64[0] >> s);
  return (int64_t)k * t0 - t1;
}

int64_t mul_kred(const uint64_t a, const uint64_t b, const uint64_t k, const size_t s){
  return kred(mul128(a, b), k, s);
}

uint64_t redc_pm(uint128_t c, const uint64_t k, const size_t s){
  uint64_t t0, t1;
  t0 = c.d64[0] & ((1UL << s) - 1UL);
  t1 = (c.d64[1] << (WORD_SIZE - s)) | (c.d64[0] >> s);
  c = mul128(t1, k);
  c = adds128_64(c, t0);
  t0 = c.d64[0] & ((1UL << s) - 1UL);
  t1 = (c.d64[1] << (WORD_SIZE - s)) | (c.d64[0] >> s);
  t0 += t1 * k;
#ifndef REDC_LAZY
  t0 = (t0 >= p_pm) ? t0 - p_pm : t0;
#endif
  return t0;
}

uint64_t redc_pm2(uint128_t c, const uint64_t k, const size_t s){
  uint64_t q, r;
  q = (c.d64[1] << (WORD_SIZE - s)) | (c.d64[0] >> s);
  r = c.d64[0] & ((1UL << s) - 1UL);
  c = mul128(q, k);
  q = (c.d64[1] << (WORD_SIZE - s)) | (c.d64[0] >> s);
  r += c.d64[0] & ((1UL << s) - 1UL);
  r += q * k;
#ifndef REDC_LAZY
  r = (r >= p_pm) ? r - p_pm : r;
#endif
  return r;
}

uint64_t mul_pm(const uint64_t a, const uint64_t b, const uint64_t k, const uint64_t s){
  return redc_pm(mul128(a, b), k, s);
}

uint64_t mul_pm2(const uint64_t a, const uint64_t b, const uint64_t k, const uint64_t s){
  return redc_pm2(mul128(a, b), k, s);
}

uint64_t redc_nfl(uint128_t c, const uint64_t p, const uint64_t v, const size_t s){
  uint128_t t;
  uint64_t t0, t1, r;
  t0 = c.d64[0];
  t1 = c.d64[1];
  c = lshift128(c, s);
  t = mul128(v, t1);
  t = adds128(t, c);
  r = t0 - t.d64[1] * p;
#ifndef REDC_LAZY
  r = (r >= p) ? r - p : r;
#endif
  return r;
}


uint64_t redc_brt(uint128_t c, const uint64_t p, const uint64_t *mu){
  uint128_t t, t_;
  // supposing that mu[0] and mu[1] are small (<= 62-bit) then there's no carry during additions.
  t = mul128(c.d64[0], mu[0]);
  t_ = mul128(c.d64[0], mu[1]);
  t_ = adds128_64(t_, t.d64[1]);
  t = mul128(c.d64[1], mu[0]);
  t_ = adds128(t_, t);
  t = mul128(c.d64[1], mu[1]);
  t = adds128_64(t, t_.d64[1]);

  t_ = mul128(p, t.d64[0]);
  t.d64[0] = 0UL;
  t.d64[1] *= p;
  t = adds128(t, t_);

  c = subs128(c, t);
#ifndef REDC_LAZY
  c.d64[0] = (c.d64[0] >= p) ? c.d64[0] - p : c.d64[0];
#endif
  return c.d64[0];
}

uint64_t redc_brt2(uint128_t c, const uint64_t p, const uint64_t mu){
  uint128_t t, t_;

  t_ = mul128(c.d64[0], mu);
  t = mul128(c.d64[1], mu);
  t = adds128_64(t, t_.d64[1]);
  t = rshift128(t, 2 * s_brt - WORD_SIZE);

  t_ = mul128(p, t.d64[0]);
  t.d64[0] = 0UL;
  t.d64[1] *= p;
  t = adds128(t, t_);
  t.d64[1] &= (1UL << (2 * s_brt - WORD_SIZE)) - 1UL;

  c = subs128(c, t);
#ifndef REDC_LAZY
  c.d64[0] = (c.d64[0] >= p) ? c.d64[0] - p : c.d64[0];
  c.d64[0] = (c.d64[0] >= p) ? c.d64[0] - p : c.d64[0];
#endif
  return c.d64[0];
}

uint64_t mul_brt(const uint64_t a, const uint64_t b, const uint64_t p, const uint64_t *mu){
  return redc_brt(mul128(a, b), p, mu);
}

uint64_t mul_brt2(const uint64_t a, const uint64_t b, const uint64_t p, const uint64_t mu){
  return redc_brt2(mul128(a, b), p, mu);
}

uint64_t mul_brt_prep(const uint64_t b, const uint64_t p){
  return ((unsigned __int128)b << 64U) / p;
}

inline uint64_t mul_brt_fixed(const uint64_t a, const uint64_t b, const uint64_t bp, const uint64_t p){
  uint128_t t;
  uint64_t r;
  t = mul128(a, bp);
  r = a * b - t.d64[1] * p;
#ifndef REDC_LAZY
  r = (r >= p) ? r - p : r;
#endif
  return r;
}

uint64_t mul_nfl(const uint64_t a, const uint64_t b, const uint64_t p, const uint64_t v, const size_t s){
  return redc_nfl(mul128(a, b), p, v, s);
}

uint64_t mul_mm(const uint64_t a, const uint64_t b, const uint64_t p, const uint64_t m){
  uint128_t t, c;
  uint64_t u, r;
  t = mul128(a, b);
  u = t.d64[0] * m;
  c = mul128(u, p);
  c = adds128(c, t);
  r = c.d64[1];
#ifndef REDC_LAZY
  r = (r >= p) ? (r - p) : r;
#endif
  return r;
}

#endif //REDUCE_GCC_H
