//
// Created by Takuya HAYASHI on 2019/11/16.
//

#ifndef FASTERNTT_FNTT_AVX2_HPP
#define FASTERNTT_FNTT_AVX2_HPP

#include "const.hpp"

#include "arith.hpp"
#include "arith_asm.hpp"
#include "arith_avx2.hpp"
#include "fntt.hpp"

#include <immintrin.h>

namespace fntt {
  namespace core {
    namespace FNTT {

      inline void forward_butterfly4_lazy_avx2(uint64_t *x, uint64_t *y, const __m256i w, const __m256i wp,
                                               const __m256i Lp, const __m256i p) {
        __m256i u, v, x_, y_;
        u = _mm256_load_si256((__m256i *) x);
        v = _mm256_load_si256((__m256i *) y);
        v = _mm256_mul_barrett_fixed_lazy(v, w, wp, p);
        x_ = _mm256_add_epi64(u, v);
        y_ = _mm256_add_epi64(u, Lp);
        y_ = _mm256_sub_epi64(y_, v);
        _mm256_store_si256((__m256i *) x, x_);
        _mm256_store_si256((__m256i *) y, y_);
      }

      inline void forward_butterfly4_lazy_avx2_2nd_to_last(uint64_t *x, uint64_t *y, const __m256i w, const __m256i wp,
                                                           const __m256i Lp, const __m256i p) {
        __m256i u, v, x_, y_;
        x_ = _mm256_load_si256((__m256i *) x); // (y1, y0, x1, x0)
        y_ = _mm256_load_si256((__m256i *) y); // (y3, y2, x3, x2)
        u = _mm256_permute4x64_epi64(y_, shuffle_imm(0, 1, 0, 1)); // (x3, x2, x3, x2)
        v = _mm256_permute4x64_epi64(x_, shuffle_imm(2, 3, 2, 3)); // (y1, y0, y1, y0)
        u = _mm256_blend_epi32(u, x_, blend32_imm(1, 1, 1, 1, 0, 0, 0, 0)); // (x3, x2, x1, x0)
        v = _mm256_blend_epi32(v, y_, blend32_imm(0, 0, 0, 0, 1, 1, 1, 1)); // (y3, y2, y1, y0)

        v = _mm256_mul_barrett_fixed_lazy(v, w, wp, p);
        x_ = _mm256_add_epi64(u, v);
        y_ = _mm256_add_epi64(u, Lp);
        y_ = _mm256_sub_epi64(y_, v);

        u = _mm256_permute4x64_epi64(x_, shuffle_imm(2, 3, 2, 3)); // (x3, x2, x3, x2)
        v = _mm256_permute4x64_epi64(y_, shuffle_imm(0, 1, 0, 1)); // (y1, y0, y1, y0)
        x_ = _mm256_blend_epi32(x_, v, blend32_imm(0, 0, 0, 0, 1, 1, 1, 1)); // (y1, y0, x1, x0)
        y_ = _mm256_blend_epi32(y_, u, blend32_imm(1, 1, 1, 1, 0, 0, 0, 0)); // (y3, y2, x3, x2)

        _mm256_store_si256((__m256i *) x, x_);
        _mm256_store_si256((__m256i *) y, y_);
      }

      inline void forward_butterfly4_montgomery_domain_lazy_avx2_last(uint64_t *x, uint64_t *y,
                                                                      const __m256i w, const __m256i wp,
                                                                      const __m256i Lp, const __m256i r,
                                                                      const __m256i rp,
                                                                      const __m256i p) {
        __m256i u, v, x_, y_;
        x_ = _mm256_load_si256((__m256i *) x); // y1 x1 y0 x0
        y_ = _mm256_load_si256((__m256i *) y); // y3 x3 y2 x2
        u = _mm256_permute4x64_epi64(x_, shuffle_imm(0, 2, 0, 2)); // x1 x0 x1 x0
        v = _mm256_permute4x64_epi64(y_, shuffle_imm(0, 2, 0, 2)); // x3 x2 x3 x2
        u = _mm256_blend_epi32(u, v, blend32_imm(0, 0, 0, 0, 1, 1, 1, 1)); // x3 x2 x1 x0
        v = _mm256_permute4x64_epi64(x_, shuffle_imm(1, 3, 1, 3)); // y1 y0 y1 y0
        y_ = _mm256_permute4x64_epi64(y_, shuffle_imm(1, 3, 1, 3)); // y3 y2 y3 y2
        v = _mm256_blend_epi32(v, y_, blend32_imm(0, 0, 0, 0, 1, 1, 1, 1)); // y3 y2 y1 y0

        u = _mm256_mul_barrett_fixed_lazy(u, r, rp, p);
        v = _mm256_mul_barrett_fixed_lazy(v, w, wp, p);
        x_ = _mm256_add_epi64(u, v);
        y_ = _mm256_add_epi64(u, Lp);
        y_ = _mm256_sub_epi64(y_, v);

        u = _mm256_permute4x64_epi64(y_, shuffle_imm(0, 0, 1, 1)); // y1 y1 y0 y0
        v = _mm256_permute4x64_epi64(x_, shuffle_imm(2, 2, 3, 3)); // x3 x3 x2 x2
        x_ = _mm256_permute4x64_epi64(x_, shuffle_imm(0, 0, 1, 1)); // x1 x1 x0 x0
        x_ = _mm256_blend_epi32(x_, u, blend32_imm(0, 0, 1, 1, 0, 0, 1, 1)); // y1 x1 y0 x0
        y_ = _mm256_permute4x64_epi64(y_, shuffle_imm(2, 2, 3, 3)); // y3 y3 y2 y2
        y_ = _mm256_blend_epi32(y_, v, blend32_imm(1, 1, 0, 0, 1, 1, 0, 0)); // y3 x3 y2 x2
        _mm256_store_epi64((__m256i *) x, x_);
        _mm256_store_epi64((__m256i *) y, y_);
      }

      inline void backward_butterfly4_lazy(uint64_t *x, uint64_t *y, const __m256i w, const __m256i wp,
                                           const __m256i Lp, const __m256i p) {
        __m256i u, v, x_, y_;
        u = _mm256_load_si256((__m256i *) x);
        v = _mm256_load_si256((__m256i *) y);
        x_ = _mm256_add_epi64(u, v);
        y_ = _mm256_add_epi64(u, Lp);
        y_ = _mm256_sub_epi64(y_, v);
        y_ = _mm256_mul_barrett_fixed_lazy(y_, w, wp, p);
        _mm256_store_si256((__m256i *) x, x_);
        _mm256_store_si256((__m256i *) y, y_);
      }

      inline void backward_butterfly4_corr_lazy(uint64_t *x, uint64_t *y, const __m256i w, const __m256i wp,
                                                const __m256i Lp, const __m256i r, const __m256i p) {
        __m256i u, v, x_, y_;
        u = _mm256_load_si256((__m256i *) x);
        v = _mm256_load_si256((__m256i *) y);
        x_ = _mm256_add_epi64(u, v);
        x_ = _mm256_redc_barrett_lazy(x_, r, p);
        y_ = _mm256_add_epi64(u, Lp);
        y_ = _mm256_sub_epi64(y_, v);
        y_ = _mm256_mul_barrett_fixed_lazy(y_, w, wp, p);
        _mm256_store_si256((__m256i *) x, x_);
        _mm256_store_si256((__m256i *) y, y_);
      }

      inline void backward_butterfly4_montgomery_domain(uint64_t *x, uint64_t *y, const __m256i w, const __m256i wp,
                                                        const __m256i ni, const __m256i nip, const __m256i Lp,
                                                        const __m256i p) {
        __m256i u, v, x_, y_;
        u = _mm256_load_si256((__m256i *) x);
        v = _mm256_load_si256((__m256i *) y);
        x_ = _mm256_add_epi64(u, v);
        x_ = _mm256_mul_barrett_fixed(x_, ni, nip, p);
        y_ = _mm256_add_epi64(u, Lp);
        y_ = _mm256_sub_epi64(y_, v);
        y_ = _mm256_mul_barrett_fixed(y_, w, wp, p);
        _mm256_store_si256((__m256i *) x, x_);
        _mm256_store_si256((__m256i *) y, y_);
      }

      // AVX2 is faster on loops for t >= 4, but slower for t = 1, 2.
      template<size_t N, enum ImplSelector Impl, std::enable_if_t<Impl == AVX2, std::nullptr_t> = nullptr>
      void ntt(uint64_t *x, const uint64_t *w, const uint64_t *wp,
               const uint64_t r, const uint64_t rp, const uint64_t p) {

        const uint64_t Lp = 2 * p;
        size_t t = N / 2;

        ++w, ++wp;

        const __m256i Lpm = _mm256_set1_epi64x(Lp), pm = _mm256_set1_epi64x(p);

        for(size_t m = 1; m < N / 4; m *= 2) {
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + t;
          for(size_t i = 0; i < m; ++i) {
            const __m256i wm = _mm256_set1_epi64x(*w), wpm = _mm256_set1_epi64x(*wp);
            for(size_t j = 0; j < t; j += 4) {
              forward_butterfly4_lazy_avx2(xp0, xp1, wm, wpm, Lpm, pm);
              xp0 += 4, xp1 += 4;
            }
            xp0 += t, xp1 += t;
            ++w, ++wp;
          }
          t /= 2;
        }

        { // t = 2
          const size_t m = N / 4;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + 2;

          for(size_t i = 0; i < m; ++i) { // unrolling
            forward_butterfly_lazy<x86_64>(xp0++, xp1++, *w, *wp, Lp, p);
            forward_butterfly_lazy<x86_64>(xp0, xp1, *w, *wp, Lp, p);

            xp0 += 2 + 1, xp1 += 2 + 1; // combine the last xp0++, xp1++ to this increment
            ++w, ++wp;
          }
        }

        { // last loop, t = 1
          const size_t m = N / 2;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + 1;
          for(size_t i = 0; i < m; ++i) {
            forward_butterfly_montgomery_domain_lazy<x86_64>(xp0, xp1, *w, *wp, Lp, r, rp, p);

            xp0 += 2, xp1 += 2; // combine the last xp0++, xp1++ to this increment
            ++w, ++wp;
          }
        }

//      const __m256i rm = _mm256_set1_epi64x(r), rpm = _mm256_set1_epi64x(rp);
//      { // t = 2, AVX2 ver. slightly slower because of complicated access patterns.
//        const size_t m = N / 4 / 2;
//        uint64_t *xp0 = x;
//        uint64_t *xp1 = x + 4;
//
//        for(size_t i = 0; i < m; ++i) { // merge two loops
//          const __m256i wm = _mm256_set_epi64x(w[1], w[1], w[0], w[0]),
//            wpm = _mm256_set_epi64x(wp[1], wp[1], wp[0], wp[0]);
//
//          forward_butterfly4_lazy_avx2_2nd_to_last(xp0, xp1, wm, wpm, Lpm, pm);
//
//          xp0 += 8, xp1 += 8;
//          w += 2, wp += 2;
//        }
//      }
//      { // last loop, t = 1, AVX2 ver. slightly slower
//        const size_t m = N / 2 / 4;
//        uint64_t *xp0 = x;
//        uint64_t *xp1 = x + 4;
//        for(size_t i = 0; i < m; ++i) {
//          const __m256i wm = _mm256_load_si256((__m256i *) w), wpm = _mm256_load_si256((__m256i *) wp);
//
//          forward_butterfly4_montgomery_domain_lazy_avx2_last(xp0, xp1, wm, wpm, Lpm, rm, rpm, pm);
//
//          xp0 += 8, xp1 += 8;
//          w += 4, wp += 4;
//        }
//      }

      }


      template<size_t N, enum ImplSelector Impl, std::enable_if_t<Impl == AVX2, std::nullptr_t> = nullptr>
      void intt(uint64_t *x, const uint64_t *w, const uint64_t *wp,
                const uint64_t ni, const uint64_t nip, const uint64_t r, const uint64_t p) {

        constexpr size_t s = kWordSize - prime59::kPrimeBitSize - 1;
        constexpr size_t L = 1UL << s;
        const uint64_t Lp = L * p;
        uint64_t *xp0, *xp1;

        ++w, ++wp;

        const __m256i Lpm = _mm256_set1_epi64x(Lp), pm = _mm256_set1_epi64x(p), rm = _mm256_set1_epi64x(r),
          nim = _mm256_set1_epi64x(ni), nipm = _mm256_set1_epi64x(nip);

        { // first loop, t = 1
          constexpr size_t m = N / 2;
          xp0 = x, xp1 = x + 1;
          for (size_t i = 0; i < m; ++i) {
            // reduce u + v since input can be in [0, Lq)
            backward_butterfly_corr_lazy<Impl>(xp0, xp1, *w, *wp, Lp, r, p);
            xp0 += 2, xp1 += 2;
            ++w, ++wp;
          }
        } // now all elements are in [0, 2q)

        { // second loop, t = 2
          constexpr size_t m = N / 4;
          xp0 = x, xp1 = x + 2;
          for (size_t i = 0; i < m; ++i) {
            backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
            backward_butterfly_lazy<Impl>(xp0, xp1, *w, *wp, Lp, p);
            xp0 += 2 + 1, xp1 += 2 + 1;
            ++w, ++wp;
          }
        } // some are in [0, 4q), others are in [0, 2q)

        constexpr int logL = s;
        constexpr int logN = (int) utils::bitLength(N) - 2; // omit the 2nd and the last layers
        constexpr int Q = logN / logL, R = logN - Q * logL;

        {
          size_t t = 4;
          size_t m = N / 8;

          for (int ll = 1; ll < logL - 1; ++ll) { // there's no need of additional reductions
            xp0 = x, xp1 = x + t;
            for (size_t i = 0; i < m; ++i) {
              const __m256i wm = _mm256_set1_epi64x(*w), wpm = _mm256_set1_epi64x(*wp);
              for (size_t j = 0; j < t; j += 4) {
                backward_butterfly4_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 4, xp1 += 4;
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          xp0 = x, xp1 = x + t;
          for (size_t i = 0; i < m; ++i) {
            const __m256i wm = _mm256_set1_epi64x(*w), wpm = _mm256_set1_epi64x(*wp);
            size_t j = 0;
            for (; j < 2 * t / L; j += 4) {
              backward_butterfly4_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
              xp0 += 4, xp1 += 4;
            }
            for (; j < t; j += 4) {
              backward_butterfly4_lazy(xp0, xp1, wm, wpm, Lpm, pm);
              xp0 += 4, xp1 += 4;
            }
            xp0 += t, xp1 += t;
            ++w, ++wp;
          }
          t *= 2, m /= 2;

          for (int qq = 1; qq < Q; ++qq) {
            for (int ll = 0; ll < logL - 1; ++ll) {
              xp0 = x, xp1 = x + t;
              for (size_t i = 0; i < m; ++i) {
                const __m256i wm = _mm256_set1_epi64x(*w), wpm = _mm256_set1_epi64x(*wp);
                size_t j = 0;
                for (; j < 1 * t / L; j += 4) {
                  backward_butterfly4_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                  xp0 += 4, xp1 += 4;
                }
                for (; j < 2 * t / L; j += 4) {
                  backward_butterfly4_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                  xp0 += 4, xp1 += 4;
                }
                for (; j < t; j += 4) {
                  backward_butterfly4_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                  xp0 += 4, xp1 += 4;
                }
                xp0 += t, xp1 += t;
                ++w, ++wp;
              }
              t *= 2, m /= 2;
            }

            xp0 = x, xp1 = x + t;
            for (size_t i = 0; i < m; ++i) {
              const __m256i wm = _mm256_set1_epi64x(*w), wpm = _mm256_set1_epi64x(*wp);
              size_t j = 0;
              for (; j < 2 * t / L; j += 4) {
                backward_butterfly4_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                xp0 += 4, xp1 += 4;
              }
              for (; j < t; j += 4) {
                backward_butterfly4_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 4, xp1 += 4;
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          for (int rr = 0; rr < R; ++rr) {
            xp0 = x, xp1 = x + t;
            for (size_t i = 0; i < m; ++i) {
              const __m256i wm = _mm256_set1_epi64x(*w), wpm = _mm256_set1_epi64x(*wp);
              size_t j = 0;
              for (; j < 1 * t / L; j += 4) {
                backward_butterfly4_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 4, xp1 += 4;
              }
              for (; j < 2 * t / L; j += 4) {
                backward_butterfly4_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                xp0 += 4, xp1 += 4;
              }
              for (; j < t; j += 4) {
                backward_butterfly4_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 4, xp1 += 4;
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          { // last loop t = n / 2
            constexpr size_t m = 1;
            constexpr size_t t = N / 2;
            xp0 = x, xp1 = x + t;
            const __m256i wm = _mm256_set1_epi64x(*w), wpm = _mm256_set1_epi64x(*wp);
            for (size_t j = 0; j < t; j += 4) { // unrolling
              backward_butterfly4_montgomery_domain(xp0, xp1, wm, wpm, nim, nipm, Lpm, pm);
              xp0 += 4, xp1 += 4;
            }
          }
        }
      }

    } // namespace FNTT
  } // namespace core
} // namespace fntt

#endif //FASTERNTT_FNTT_AVX2_HPP
