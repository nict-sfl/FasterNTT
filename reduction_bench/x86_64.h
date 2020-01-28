//
// Created by Takuya HAYASHI on 2019/11/11.
//

#ifndef REDUCE_X86_64_H
#define REDUCE_X86_64_H

#include "param.h"

uint64_t mul_nfl_asm(const uint64_t a, const uint64_t b, const uint64_t p, const uint64_t v, const uint8_t s) {
  uint64_t c0, c1, t0, t1, r;
  __asm__ volatile (
  "mulx %[b], %[c0], %%rdx\n\t"
  "mov %[c0], %[r]\n\t"
  "mulx %[v], %[t0], %[t1]\n\t"
  "shld %[s], %[c0], %%rdx\n\t"
  "shl %[s], %[c0]\n\t"
  "add %[c0], %[t0]\n\t"
  "adc %%rdx, %[t1]\n\t"
  "imul %[p], %[t1]\n\t"
  "sub %[t1], %[r]\n\t"
#ifndef REDC_LAZY
  "mov %[r], %[c0]\n\t"
  "sub %[p], %[c0]\n\t"
  "cmovnc %[c0], %[r]\n\t"
#endif
  : [r] "=&r"(r), [t0] "=&r"(t0), [t1] "=&r"(t1), [c0] "=&r"(c0), [c1] "=&r"(c1)
  : [a] "d"(a), [b] "r"(b), [s]"c"(s), [v]"r"(v), [p]"r"(p)
  : "cc");
  return r;
}

int64_t mul_kred_asm(const uint64_t a, const uint64_t b, const uint64_t k, const uint8_t s) {
  const uint64_t mask = (1UL << s) - 1UL;
  uint64_t c0, c1;
  int64_t t;
  __asm__ volatile (
  "mulx %[b], %[c0], %[c1]\n\t"
  "mov %[c0], %[t]\n\t"
  "shrd %[s], %[c1], %[c0]\n\t"
  "and %[mask], %[t]\n\t"
  "imul %[k], %[t]\n\t"
  "sub %[c0], %[t]\n\t"
  : [t] "=&r"(t), [c0] "=&r"(c0), [c1] "=&r"(c1)
  : [a] "d"(a), [b] "r"(b), [mask]"r"(mask), [s]"c"(s), [k]"r"(k)
  : "cc");
  return t;
}

uint64_t mul_mm_asm(const uint64_t a, const uint64_t b, const uint64_t p, const uint64_t m) {
  uint64_t t0, t1, c0, c1;
  __asm__ volatile (
  "mulx %[b], %[t0], %[t1]\n\t"
  "mov %[m], %[c0]\n\t"
  "mov %[p], %%rdx\n\t"
  "imul %[t0], %[c0]\n\t"
  "mulx %[c0], %[c0], %[c1]\n\t"
  "add %[t0], %[c0]\n\t"
  "adc %[t1], %[c1]\n\t"
#ifndef REDC_LAZY
  "mov %[c1], %[t0]\n\t"
  "sub %[p], %[t0]\n\t"
  "cmovnc %[t0], %[c1]\n\t"
#endif
  : [c1] "=&r"(c1), [c0] "=&r"(c0), [t1] "=&r"(t1), [t0] "=&r"(t0)
  : [a] "d"(a), [b] "r"(b), [p] "r"(p), [m] "r"(m)
  : "cc");
  return c1;
}

uint64_t mul_brt_fixed_asm(uint64_t a, const uint64_t b, const uint64_t bp, const uint64_t p) {
  uint64_t t0, t1;
  __asm__ volatile (
  "mulx %[bp], %[t0], %[t1]\n\t"
  "imul %[b], %[a]\n\t"
  "imul %[p], %[t1]\n\t"
  "sub %[t1], %[a]\n\t"
#ifndef REDC_LAZY
  "mov %[a], %[t0]\n\t"
  "sub %[p], %[t0]\n\t"
  "cmovnc %[t0], %[a]\n\t"
#endif
  : [a] "+&d"(a), [t0] "=&r"(t0), [t1] "=&r"(t1)
  : [b] "r"(b), [p] "r"(p), [bp] "r"(bp)
  : "cc");
  return a;
}

#endif //REDUCE_X86_64_H
