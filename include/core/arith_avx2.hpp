//
// Created by Takuya HAYASHI on 2019/11/13.
//

#ifndef REDUCE_AVX_H
#define REDUCE_AVX_H

#include "utils_simd.hpp"
#include <immintrin.h>

namespace fntt {
  namespace core {
    inline void _mm256_add_epu128(__m256i *c1, __m256i *c0, const __m256i a1, const __m256i a0, const __m256i b1,
                                  const __m256i b0) {
      const __m256i m = _mm256_set1_epi64x(0xffffffffUL);
      __m256i s, t, u;
      s = _mm256_and_si256(a0, m);
      t = _mm256_and_si256(b0, m);
      u = _mm256_add_epi64(s, t);
      s = _mm256_srli_epi64(a0, 32);
      t = _mm256_srli_epi64(b0, 32);
      u = _mm256_add_epi64(u, s);
      u = _mm256_add_epi64(u, t);
      *c0 = _mm256_add_epi64(a0, b0);
      u = _mm256_srli_epi64(u, 32);
      *c1 = _mm256_add_epi64(a1, b1);
      *c1 = _mm256_add_epi64(*c1, u);
    }

    // a = (a >= p) ? a - p : a
    inline __m256i _mm256_conditional_sub_epu64(const __m256i a, const __m256i p) {
      const __m256i z = _mm256_setzero_si256();
      __m256i t, m;
      t = _mm256_sub_epi64(a, p);
      m = _mm256_srli_epi64(t, 63);
      m = _mm256_sub_epi64(z, m);
      m = _mm256_and_si256(m, p);
      t = _mm256_add_epi64(t, m);
      return t;
    }

    inline void _mm256_mul_epu64(__m256i *c1, __m256i *c0, const __m256i a, const __m256i b) {
      const __m256i mask = _mm256_set1_epi64x(0xffffffffUL);
      __m256i a1, b1, s, s0, t, t0, u;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm256_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      *c0 = _mm256_mul_epu32(a, b); // a0 x b0
      t = _mm256_mul_epu32(a1, b); // a1 x b0 (2^32)
      s = _mm256_mul_epu32(b1, a); // a0 x b1 (2^32)
      *c1 = _mm256_mul_epu32(a1, b1); // a1 x b1 (2^64)

      u = _mm256_srli_epi64(*c0, 32);
      t0 = _mm256_and_si256(t, mask);
      s0 = _mm256_and_si256(s, mask);
      t = _mm256_srli_epi64(t, 32);
      s = _mm256_srli_epi64(s, 32);

      u = _mm256_add_epi64(u, t0);
      u = _mm256_add_epi64(u, s0);
      t0 = _mm256_shuffle_epi32(u, shuffle_imm(0, 0, 2, 2));
      *c0 = _mm256_blend_epi32(*c0, t0, blend32_imm(0, 1, 0, 1, 0, 1, 0, 1));

      u = _mm256_srli_epi64(u, 32);
      *c1 = _mm256_add_epi64(*c1, u);
      *c1 = _mm256_add_epi64(*c1, t);
      *c1 = _mm256_add_epi64(*c1, s);
    }

    // return upper 64-bit of (a * b) + (2^64 c1 + c0)
    inline __m256i _mm256_addmulhi_epu64(const __m256i a, const __m256i b, const __m256i c1, const __m256i c0) {
      const __m256i mask = _mm256_set1_epi64x(0xffffffffUL);
      __m256i a1, b1, s, s0, t, t0, u0, u00, u1, c01, r;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm256_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      u0 = _mm256_mul_epu32(a, b); // a0 x b0
      t = _mm256_mul_epu32(a1, b); // a1 x b0 (2^32)
      s = _mm256_mul_epu32(b1, a); // a0 x b1 (2^32)
      u1 = _mm256_mul_epu32(a1, b1); // a1 x b1 (2^64)

      u00 = _mm256_and_si256(u0, mask);
      r = _mm256_and_si256(c0, mask);
      u0 = _mm256_srli_epi64(u0, 32);
      c01 = _mm256_srli_epi64(c0, 32);
      u0 = _mm256_add_epi64(u0, c01);
      r = _mm256_add_epi64(r, u00);
      r = _mm256_srli_epi64(r, 32);

      t0 = _mm256_and_si256(t, mask);
      s0 = _mm256_and_si256(s, mask);
      r = _mm256_add_epi64(r, u0);
      r = _mm256_add_epi64(r, t0);
      r = _mm256_add_epi64(r, s0);

      t = _mm256_srli_epi64(t, 32);
      s = _mm256_srli_epi64(s, 32);
      r = _mm256_srli_epi64(r, 32);

      t = _mm256_add_epi64(t, s);
      r = _mm256_add_epi64(r, u1);
      r = _mm256_add_epi64(r, t);
      r = _mm256_add_epi64(r, c1);

      return r;
    }

    inline __m256i _mm256_mulhi_epu64(const __m256i a, const __m256i b) {
      const __m256i mask = _mm256_set1_epi64x(0xffffffffUL);
      __m256i a1, b1, s, s0, t, t0, c0, c1;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm256_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      c0 = _mm256_mul_epu32(a, b); // a0 x b0
      t = _mm256_mul_epu32(a1, b); // a1 x b0 (2^32)
      s = _mm256_mul_epu32(b1, a); // a0 x b1 (2^32)
      c1 = _mm256_mul_epu32(a1, b1); // a1 x b1 (2^64)

      c0 = _mm256_srli_epi64(c0, 32);
      t0 = _mm256_and_si256(t, mask);
      s0 = _mm256_and_si256(s, mask);
      t = _mm256_srli_epi64(t, 32);
      s = _mm256_srli_epi64(s, 32);

      c0 = _mm256_add_epi64(c0, t0);
      c0 = _mm256_add_epi64(c0, s0);
      c0 = _mm256_srli_epi64(c0, 32);
      c1 = _mm256_add_epi64(c1, c0);
      c1 = _mm256_add_epi64(c1, t);
      c1 = _mm256_add_epi64(c1, s);
      return c1;
    }

    inline __m256i _mm256_mullo_epu64(const __m256i a, const __m256i b) {
      __m256i a1, b1, t, s, c0;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm256_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      c0 = _mm256_mul_epu32(a, b); // a0 x b0
      t = _mm256_mul_epu32(a1, b); // a1 x b0 (2^32)
      s = _mm256_mul_epu32(b1, a); // a0 x b1 (2^32)
      t = _mm256_add_epi64(t, s); // a1 x b0 + a0 x b1 (2^32)
      t = _mm256_slli_epi64(t, 32);
      c0 = _mm256_add_epi64(c0, t);
      return c0;
    }

    // suppose a, b < 2^62
    inline void _mm256_mul_epu64_ab62(__m256i *c1, __m256i *c0, const __m256i a, const __m256i b) {
      __m256i a1, b1, s, t, u;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm256_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      *c0 = _mm256_mul_epu32(a, b); // a0 x b0
      t = _mm256_mul_epu32(a1, b); // a1 x b0 (2^32)
      s = _mm256_mul_epu32(b1, a); // a0 x b1 (2^32)
      *c1 = _mm256_mul_epu32(a1, b1); // a1 x b1 (2^64)

      u = _mm256_srli_epi64(*c0, 32);
      t = _mm256_add_epi64(t, s);
      t = _mm256_add_epi64(t, u);
      s = _mm256_slli_epi64(t, 32);
      t = _mm256_srli_epi64(t, 32);

      *c0 = _mm256_blend_epi32(*c0, s, blend32_imm(0, 1, 0, 1, 0, 1, 0, 1));
      *c1 = _mm256_add_epi64(*c1, t);
    }

    // suppose b < 2^62
    inline __m256i _mm256_mulhi_epu64_b62(const __m256i a, const __m256i b) {
      const __m256i mask = _mm256_set1_epi64x(0xffffffffUL);
      __m256i a1, b1, s, s0, t, c0, c1;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm256_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      c0 = _mm256_mul_epu32(a, b); // a0b0
      t = _mm256_mul_epu32(a1, b); // a1b0
      s = _mm256_mul_epu32(b1, a); // a0b1
      c1 = _mm256_mul_epu32(a1, b1); // a1b1

      c0 = _mm256_srli_epi64(c0, 32);
      s0 = _mm256_and_si256(s, mask);
      s = _mm256_srli_epi64(s, 32);
      t = _mm256_add_epi64(t, c0);
      t = _mm256_add_epi64(t, s0);
      t = _mm256_srli_epi64(t, 32);

      c1 = _mm256_add_epi64(c1, s);
      c1 = _mm256_add_epi64(c1, t);
      return c1;
    }

    // suppose b < 2^59
    inline __m256i _mm256_addmulhi_epu64_b59(const __m256i a, const __m256i b, const __m256i c1, const __m256i c0) {
      const __m256i mask = _mm256_set1_epi64x(0xffffffffUL);
      __m256i a1, b1, s, s0, t, t0, u0, u00, u1, c01, r;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm256_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      u0 = _mm256_mul_epu32(a, b); // a0 x b0
      t = _mm256_mul_epu32(a1, b); // a1 x b0 (2^32)
      s = _mm256_mul_epu32(b1, a); // a0 x b1 (2^32)
      u1 = _mm256_mul_epu32(a1, b1); // a1 x b1 (2^64)

      u00 = _mm256_and_si256(u0, mask);
      r = _mm256_and_si256(c0, mask);
      u0 = _mm256_srli_epi64(u0, 32);
      c01 = _mm256_srli_epi64(c0, 32);
      r = _mm256_add_epi64(r, u00);
      u0 = _mm256_add_epi64(u0, c01);
      r = _mm256_srli_epi64(r, 32);
      r = _mm256_add_epi64(r, u0);
      r = _mm256_add_epi64(r, s);
      t0 = _mm256_and_si256(t, mask);
      r = _mm256_add_epi64(r, t0);
      t = _mm256_srli_epi64(t, 32);
      r = _mm256_srli_epi64(r, 32);
      r = _mm256_add_epi64(r, t);
      r = _mm256_add_epi64(r, c1);
      r = _mm256_add_epi64(r, u1);
      return r;
    }

    // suppose b < 2^32
    inline __m256i _mm256_mulhi_epu64_b32(const __m256i a, const __m256i b) {
      __m256i a1, c0, c1;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      c0 = _mm256_mul_epu32(a, b);
      c1 = _mm256_mul_epu32(a1, b);
      c0 = _mm256_srli_epi64(c0, 32);
      c1 = _mm256_add_epi64(c1, c0);
      c1 = _mm256_srli_epi64(c1, 32);
      return c1;
    }

    // suppose b < 2^32
    inline __m256i _mm256_mullo_epu64_b32(const __m256i a, const __m256i b) {
      __m256i a1, c0, c1;
      a1 = _mm256_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      c1 = _mm256_mul_epu32(a1, b);
      c0 = _mm256_mul_epu32(a, b);
      c1 = _mm256_slli_epi64(c1, 32);
      c0 = _mm256_add_epi64(c0, c1);
      return c0;
    }

    inline __m256i _mm256_redc_barrett_lazy(const __m256i a, const __m256i r, const __m256i p) {
      __m256i t;
      t = _mm256_mulhi_epu64_b32(a, r);
      t = _mm256_mullo_epu64_b32(p, t);
      t = _mm256_sub_epi64(a, t);
      return t;
    }

    inline __m256i _mm256_redc_barrett(const __m256i a, const __m256i r, const __m256i p) {
      __m256i t;
      t = _mm256_mulhi_epu64_b32(a, r);
      t = _mm256_mullo_epu64_b32(p, t);
      t = _mm256_sub_epi64(a, t);
      t = _mm256_conditional_sub_epu64(t, p);
      return t;
    }

    // suppose b < 2^59
    inline __m256i _mm256_mul_barrett_fixed_lazy(const __m256i a, const __m256i b, const __m256i bp,
                                                 const __m256i p) {
      __m256i t, c;
      t = _mm256_mulhi_epu64_b62(bp, a);
      c = _mm256_mullo_epu64(a, b);
      t = _mm256_mullo_epu64(t, p);
      c = _mm256_sub_epi64(c, t);
      return c;
    }

    // suppose b < 2^59
    inline __m256i _mm256_mul_barrett_fixed(const __m256i a, const __m256i b, const __m256i bp,
                                            const __m256i p) {
      __m256i t, c;
      t = _mm256_mulhi_epu64_b62(bp, a);
      c = _mm256_mullo_epu64(a, b);
      t = _mm256_mullo_epu64(t, p);
      c = _mm256_sub_epi64(c, t);
      c = _mm256_conditional_sub_epu64(c, p);
      return c;
    }

    inline __m256i _mm256_mul_montgomery_lazy(const __m256i a, const __m256i b, const __m256i m, const __m256i p) {
      __m256i t1, t0, c1, u;
      _mm256_mul_epu64_ab62(&t1, &t0, a, b);
      u = _mm256_mullo_epu64(t0, m);
      c1 = _mm256_addmulhi_epu64_b59(u, p, t1, t0);
      return c1;
    }

    inline __m256i _mm256_mul_montgomery(const __m256i a, const __m256i b, const __m256i m, const __m256i p) {
      __m256i t1, t0, c1, u;
      _mm256_mul_epu64_ab62(&t1, &t0, a, b);
      u = _mm256_mullo_epu64(t0, m);
      c1 = _mm256_addmulhi_epu64_b59(u, p, t1, t0);
      c1 = _mm256_conditional_sub_epu64(c1, p);
      return c1;
    }

  } // namespace core
} // namespace fntt

#endif //REDUCE_AVX_H
