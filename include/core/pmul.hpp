/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_PMUL_HPP
#define FASTERNTT_PMUL_HPP

#include "const.hpp"
#include "arith.hpp"
#include "arith_asm.hpp"

namespace fntt {
  namespace core {

    /// Input should be in Montgomery domain, and a[i] * b[i] must be in [0, Rq).
    /// Output will be in [0, q)
    template<size_t N, enum ImplSelector Impl>
    static void pmul_montgomery(uint64_t *c, const uint64_t *a, const uint64_t *b,
                                const uint64_t m, const uint64_t p) {
      for(size_t i = 0; i < N; ++i) {
        c[i] = mul_montgomery<Impl>(a[i], b[i], m, p);
      }
    }

    /// Input should be in Montgomery domain, and a[i] * b[i] must be in [0, Rq).
    /// Output will be in [0, 2q)
    template<size_t N, enum ImplSelector Impl>
    static void pmul_montgomery_lazy(uint64_t *c, const uint64_t *a, const uint64_t *b,
                                     const uint64_t m, const uint64_t p) {
      for(size_t i = 0; i < N; ++i) {
        c[i] = mul_montgomery_lazy<Impl>(a[i], b[i], m, p);
      }
    }

//  AVX2 ver. much slower...
//    template<size_t N>
//    static void pmul(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t m, const uint64_t p) {
//      const __m256i m_ = _mm256_set1_epi64x(m), p_ = _mm256_set1_epi64x(p);
//      for(size_t i = 0; i < N; i += 4) {
//        __m256i a_ = _mm256_load_si256((__m256i *) (a + i));
//        __m256i b_ = _mm256_load_si256((__m256i *) (b + i));
//        __m256i c_ = _mm256_mul_montgomery(a_, b_, m_, p_);
//        _mm256_store_si256((__m256i *) (c + i), c_);
//      }
//  }

//  AVX2 ver. much slower...
//    template<size_t N>
//    static void pmul_lazy(uint64_t *c, const uint64_t *a, const uint64_t *b, const uint64_t m, const uint64_t p) {
//      const __m256i m_ = _mm256_set1_epi64x(m), p_ = _mm256_set1_epi64x(p);
//      for(size_t i = 0; i < N; i += 4) {
//        __m256i a_ = _mm256_load_si256((__m256i *) (a + i));
//        __m256i b_ = _mm256_load_si256((__m256i *) (b + i));
//        __m256i c_ = _mm256_mul_montgomery_lazy(a_, b_, m_, p_);
//        _mm256_store_si256((__m256i *) (c + i), c_);
//      }
//  }

    template<size_t N, enum ImplSelector Impl>
    void pmul_nfl(uint64_t *c, const uint64_t *a, const uint64_t *b,
                  const uint64_t v, const size_t s, const uint64_t p) {
      for(size_t i = 0; i < N; ++i) {
        c[i] = mul_nfl<Impl>(a[i], b[i], v, s, p);
      }
    }

    template<size_t N, enum ImplSelector Impl>
    void pmul_lazy_nfl(uint64_t *c, const uint64_t *a, const uint64_t *b,
                       const uint64_t v, const size_t s, const uint64_t p) {
      for(size_t i = 0; i < N; ++i) {
        c[i] = mul_nfl_lazy<Impl>(a[i], b[i], v, s, p);
      }
    }

  } // namespace core
} // namespace fntt

#endif //FASTERNTT_PMUL_HPP
