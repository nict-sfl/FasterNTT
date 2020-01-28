//
// Created by Takuya HAYASHI on 2019/11/11.
//

#ifndef REDUCE_PARAM_H
#define REDUCE_PARAM_H

#define WORD_SIZE 64
#define REDC_LAZY

#include <stddef.h>
#include <stdint.h>


// p = k * 2^s + 1
static const uint64_t p_kred = 180143985094819841UL;
static const uint64_t k_kred = 5;
static const uint64_t k__kred = 144115188075855873UL;
static const size_t s_kred = 55;

// p = 2^s - k, k = ell * 2^d - 1
static const uint64_t p_pm = 1152921504577486849UL;
static const uint64_t k_pm = 29360127UL;
static const size_t s_pm = 60;

// NFLlib primes
static const uint64_t p_nfl = 4611686018326724609UL; // 62-bit
static const size_t s0 = 2; // WORD_SIZE - log_2(p)
static const uint64_t v0 = 1610612720UL;

// Montgomery param
static const uint64_t p_mm = 1152921504577486849UL; // p_pm
static const uint64_t m_p_mm_i = 1152059487461310463UL; // -p_mm^-1 mod 2^WORD_SIZE
static const uint64_t rr = 220676366708769024UL; // 2^128 mod p_mm

// Barrett param
static const uint64_t p_brt = 1152921504577486849UL;
static const uint64_t mu_brt[2] = {7516192512, 16};

static const uint64_t p_brt2 = 1152921504577486849UL;
static const uint64_t mu_brt2 = 1152921504636207103UL;
static const size_t s_brt = 60;

#endif //REDUCE_PARAM_H
