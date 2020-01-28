/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_FNTT_SSE_HPP
#define FASTERNTT_FNTT_SSE_HPP

#include "const.hpp"
#include "arith.hpp"
#include "arith_asm.hpp"
#include "arith_sse.hpp"
#include "fntt.hpp"

#include <emmintrin.h>

namespace fntt {
  namespace core {
    namespace FNTT {
      inline void forward_butterfly2_lazy_sse(uint64_t *x, uint64_t *y, const __m128i w, const __m128i wp,
                                              const __m128i Lp, const __m128i p) {
        __m128i u, v, x_, y_;
        u = _mm_load_si128((__m128i *) x);
        v = _mm_load_si128((__m128i *) y);
        v = _mm_mul_barrett_fixed_lazy(v, w, wp, p);
        x_ = _mm_add_epi64(u, v);
        y_ = _mm_add_epi64(u, Lp);
        y_ = _mm_sub_epi64(y_, v);
        _mm_store_si128((__m128i *) x, x_);
        _mm_store_si128((__m128i *) y, y_);
      }

      inline void backward_butterfly2_lazy(uint64_t *x, uint64_t *y, const __m128i w, const __m128i wp,
                                           const __m128i Lp, const __m128i p) {
        __m128i u, v, x_, y_;
        u = _mm_load_si128((__m128i *) x);
        v = _mm_load_si128((__m128i *) y);
        x_ = _mm_add_epi64(u, v);
        y_ = _mm_add_epi64(u, Lp);
        y_ = _mm_sub_epi64(y_, v);
        y_ = _mm_mul_barrett_fixed_lazy(y_, w, wp, p);
        _mm_store_si128((__m128i *) x, x_);
        _mm_store_si128((__m128i *) y, y_);
      }

      inline void backward_butterfly2_corr_lazy(uint64_t *x, uint64_t *y, const __m128i w, const __m128i wp,
                                                const __m128i Lp, const __m128i r, const __m128i p) {
        __m128i u, v, x_, y_;
        u = _mm_load_si128((__m128i *) x);
        v = _mm_load_si128((__m128i *) y);
        x_ = _mm_add_epi64(u, v);
        x_ = _mm_redc_barrett_lazy(x_, r, p);
        y_ = _mm_add_epi64(u, Lp);
        y_ = _mm_sub_epi64(y_, v);
        y_ = _mm_mul_barrett_fixed_lazy(y_, w, wp, p);
        _mm_store_si128((__m128i *) x, x_);
        _mm_store_si128((__m128i *) y, y_);
      }

      inline void backward_butterfly2_montgomery_domain(uint64_t *x, uint64_t *y, const __m128i w, const __m128i wp,
                                                        const __m128i ni, const __m128i nip, const __m128i Lp,
                                                        const __m128i p) {
        __m128i u, v, x_, y_;
        u = _mm_load_si128((__m128i *) x);
        v = _mm_load_si128((__m128i *) y);
        x_ = _mm_add_epi64(u, v);
        x_ = _mm_mul_barrett_fixed(x_, ni, nip, p);
        y_ = _mm_add_epi64(u, Lp);
        y_ = _mm_sub_epi64(y_, v);
        y_ = _mm_mul_barrett_fixed(y_, w, wp, p);
        _mm_store_si128((__m128i *) x, x_);
        _mm_store_si128((__m128i *) y, y_);
      }

      template<size_t N, enum ImplSelector Impl, std::enable_if_t<Impl == SSE, std::nullptr_t> = nullptr>
      void ntt(uint64_t *x, const uint64_t *w, const uint64_t *wp,
               const uint64_t r, const uint64_t rp, const uint64_t p) {

        const uint64_t Lp = 2 * p;
        size_t t = N / 2;

        ++w, ++wp;

        const __m128i Lpm = _mm_set1_epi64x(Lp), pm = _mm_set1_epi64x(p);

        for(size_t m = 1; m < N / 4; m *= 2) {
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + t;
          for(size_t i = 0; i < m; ++i) {
            const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
            for(size_t j = 0; j < t; j += 4) {
              forward_butterfly2_lazy_sse(xp0, xp1, wm, wpm, Lpm, pm);
              xp0 += 2, xp1 += 2;
              forward_butterfly2_lazy_sse(xp0, xp1, wm, wpm, Lpm, pm);
              xp0 += 2, xp1 += 2;
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
            const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
            forward_butterfly2_lazy_sse(xp0, xp1, wm, wpm, Lpm, pm);

            xp0 += 4, xp1 += 4;
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
      }


      template<size_t N, enum ImplSelector Impl, std::enable_if_t<Impl == SSE, std::nullptr_t> = nullptr>
      void intt(uint64_t *x, const uint64_t *w, const uint64_t *wp,
                const uint64_t ni, const uint64_t nip, const uint64_t r, const uint64_t p) {

        constexpr size_t s = kWordSize - prime59::kPrimeBitSize - 1;
        constexpr size_t L = 1UL << s;
        const uint64_t Lp = L * p;
        uint64_t *xp0, *xp1;

        ++w, ++wp;

        const __m128i Lpm = _mm_set1_epi64x(Lp), pm = _mm_set1_epi64x(p), rm = _mm_set1_epi64x(r),
          nim = _mm_set1_epi64x(ni), nipm = _mm_set1_epi64x(nip);

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
              const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
              for (size_t j = 0; j < t; j += 4) {
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          xp0 = x, xp1 = x + t;
          for (size_t i = 0; i < m; ++i) {
            const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
            size_t j = 0;
            for (; j < 2 * t / L; j += 4) {
              backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
              xp0 += 2, xp1 += 2;
              backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
              xp0 += 2, xp1 += 2;
            }
            for (; j < t; j += 4) {
              backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
              xp0 += 2, xp1 += 2;
              backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
              xp0 += 2, xp1 += 2;
            }
            xp0 += t, xp1 += t;
            ++w, ++wp;
          }
          t *= 2, m /= 2;

          for (int qq = 1; qq < Q; ++qq) {
            for (int ll = 0; ll < logL - 1; ++ll) {
              xp0 = x, xp1 = x + t;
              for (size_t i = 0; i < m; ++i) {
                const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
                size_t j = 0;
                for (; j < 1 * t / L; j += 4) {
                  backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                  xp0 += 2, xp1 += 2;
                  backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                  xp0 += 2, xp1 += 2;
                }
                for (; j < 2 * t / L; j += 4) {
                  backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                  xp0 += 2, xp1 += 2;
                  backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                  xp0 += 2, xp1 += 2;
                }
                for (; j < t; j += 4) {
                  backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                  xp0 += 2, xp1 += 2;
                  backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                  xp0 += 2, xp1 += 2;
                }
                xp0 += t, xp1 += t;
                ++w, ++wp;
              }
              t *= 2, m /= 2;
            }

            xp0 = x, xp1 = x + t;
            for (size_t i = 0; i < m; ++i) {
              const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
              size_t j = 0;
              for (; j < 2 * t / L; j += 4) {
                backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                xp0 += 2, xp1 += 2;
                backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                xp0 += 2, xp1 += 2;
              }
              for (; j < t; j += 4) {
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          for (int rr = 0; rr < R; ++rr) {
            xp0 = x, xp1 = x + t;
            for (size_t i = 0; i < m; ++i) {
              const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
              size_t j = 0;
              for (; j < 1 * t / L; j += 4) {
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
              }
              for (; j < 2 * t / L; j += 4) {
                backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                xp0 += 2, xp1 += 2;
                backward_butterfly2_corr_lazy(xp0, xp1, wm, wpm, Lpm, rm, pm);
                xp0 += 2, xp1 += 2;
              }
              for (; j < t; j += 4) {
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
                backward_butterfly2_lazy(xp0, xp1, wm, wpm, Lpm, pm);
                xp0 += 2, xp1 += 2;
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          { // last loop t = n / 2
            // constexpr size_t m = 1;
            constexpr size_t t = N / 2;
            xp0 = x, xp1 = x + t;
            const __m128i wm = _mm_set1_epi64x(*w), wpm = _mm_set1_epi64x(*wp);
            for (size_t j = 0; j < t; j += 4) { // unrolling
              backward_butterfly2_montgomery_domain(xp0, xp1, wm, wpm, nim, nipm, Lpm, pm);
              xp0 += 2, xp1 += 2;
              backward_butterfly2_montgomery_domain(xp0, xp1, wm, wpm, nim, nipm, Lpm, pm);
              xp0 += 2, xp1 += 2;
            }
          }
        }
      }

    } // namespace FNTT
  } // namespace core
} // namespace fntt

#endif //FASTERNTT_FNTT_SSE_HPP
