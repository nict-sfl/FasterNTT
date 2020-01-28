/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_ARITH_SSE_HPP
#define FASTERNTT_ARITH_SSE_HPP

#include "utils_simd.hpp"
#include <emmintrin.h>

namespace fntt {
  namespace core {

    inline __m128i _mm_conditional_sub_epu64(const __m128i a, const __m128i p){
      const __m128i z = _mm_setzero_si128();
      __m128i t, m;
      t = _mm_sub_epi64(a, p);
      m = _mm_srli_epi64(t, 63);
      m = _mm_sub_epi64(z, m);
      m = _mm_and_si128(m, p);
      t = _mm_add_epi64(t, m);
      return t;
    }

    inline __m128i _mm_mulhi_epu64_b62(const __m128i a, const __m128i b){
      const __m128i mask = _mm_set1_epi64x(0xffffffffUL);
      __m128i a1, b1, s, s0, t, c0, c1;
      a1 = _mm_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      c0 = _mm_mul_epu32(a, b);
      t = _mm_mul_epu32(a1, b);
      s = _mm_mul_epu32(b1, a);
      c1 = _mm_mul_epu32(a1, b1);

      c0 = _mm_srli_epi64(c0, 32);
      s0 = _mm_and_si128(s, mask);
      s = _mm_srli_epi64(s, 32);
      t = _mm_add_epi64(t, c0);
      t = _mm_add_epi64(t, s0);
      t = _mm_srli_epi64(t, 32);

      c1 = _mm_add_epi64(c1, s);
      c1 = _mm_add_epi64(c1, t);
      return c1;
    }

    inline __m128i _mm_mullo_epu64(const __m128i a, const __m128i b){
      __m128i a1, b1, t, s, c0;
      a1 = _mm_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      b1 = _mm_shuffle_epi32(b, shuffle_imm(1, 1, 3, 3));

      c0 = _mm_mul_epu32(a, b);
      t = _mm_mul_epu32(a1, b);
      s = _mm_mul_epu32(b1, a);
      t = _mm_add_epi64(t, s);
      t = _mm_slli_epi64(t, 32);
      c0 = _mm_add_epi64(c0, t);
      return c0;
    }

    inline __m128i _mm_mulhi_epu64_b32(const __m128i a, const __m128i b){
      __m128i a1, c0, c1;
      a1 = _mm_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      c0 = _mm_mul_epu32(a, b);
      c1 = _mm_mul_epu32(a1, b);
      c0 = _mm_srli_epi64(c0, 32);
      c1 = _mm_add_epi64(c1, c0);
      c1 = _mm_srli_epi64(c1, 32);
      return c1;
    }

    inline __m128i _mm_mullo_epu64_b32(const __m128i a, const __m128i b){
      __m128i a1, c0, c1;
      a1 = _mm_shuffle_epi32(a, shuffle_imm(1, 1, 3, 3));
      c1 = _mm_mul_epu32(a1, b);
      c0 = _mm_mul_epu32(a, b);
      c1 = _mm_slli_epi64(c1, 32);
      c0 = _mm_add_epi64(c0, c1);
      return c0;
    }

    inline __m128i _mm_redc_barrett_lazy(const __m128i a, const __m128i r, const __m128i p){
      __m128i t;
      t = _mm_mulhi_epu64_b32(a, r);
      t = _mm_mullo_epu64_b32(p, t);
      t = _mm_sub_epi64(a, t);
      return t;
    }

    inline __m128i _mm_redc_barrett(const __m128i a, const __m128i r, const __m128i p){
      __m128i t;
      t = _mm_mulhi_epu64_b32(a, r);
      t = _mm_mullo_epu64_b32(p, t);
      t = _mm_sub_epi64(a, t);
      t = _mm_conditional_sub_epu64(t, p);
      return t;
    }

    inline __m128i _mm_mul_barrett_fixed_lazy(const __m128i a, const __m128i b, const __m128i bp, const __m128i p){
      __m128i t, c;
      t = _mm_mulhi_epu64_b62(bp, a);
      c = _mm_mullo_epu64(a, b);
      t = _mm_mullo_epu64(t, p);
      c = _mm_sub_epi64(c, t);
      return c;
    }

    inline __m128i _mm_mul_barrett_fixed(const __m128i a, const __m128i b, const __m128i bp, const __m128i p){
      __m128i t, c;
      t = _mm_mulhi_epu64_b62(bp, a);
      c = _mm_mullo_epu64(a, b);
      t = _mm_mullo_epu64(t, p);
      c = _mm_sub_epi64(c, t);
      c = _mm_conditional_sub_epu64(c, p);
      return c;
    }

  } // namespace core
} // namespace fntt

#endif //FASTERNTT_ARITH_SSE_HPP
