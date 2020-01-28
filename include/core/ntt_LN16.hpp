/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_NTT_LN16_HPP
#define FASTERNTT_NTT_LN16_HPP

#include "const.hpp"
#include "arith.hpp"
#include "arith_asm.hpp"

namespace fntt {
  namespace core {
    namespace LN16 {
      template <enum ImplSelector Impl>
      inline void forward_butterfly_lazy(uint64_t *x, uint64_t *y,
                                         const uint64_t w, const uint64_t wp, const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = sub_conditional<Impl>(*x, Lp);
        v = mul_barrett_fixed_lazy<Impl>(*y, w, wp, p);
        *x = u + v;
        *y = u + Lp - v;
      }

      template <enum ImplSelector Impl>
      inline void backward_butterfly_lazy(uint64_t *x, uint64_t *y,
                                          const uint64_t w, const uint64_t wp, const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = *x;
        v = *y;
        *x = u + v;
        *x = sub_conditional<Impl>(*x, Lp);
        *y = mul_barrett_fixed_lazy<Impl>(u + Lp - v, w, wp, p);
      }

      template <enum ImplSelector Impl>
      inline void backward_butterfly_last(uint64_t *x, uint64_t *y, const uint64_t w, const uint64_t wp,
                                          const uint64_t ni, const uint64_t nip, const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = *x;
        v = *y;
        *x = u + v;
        *x = sub_conditional<Impl>(*x, Lp);
        *x = mul_barrett_fixed<Impl>(*x, ni, nip, p);
        *y = mul_barrett_fixed<Impl>(u + Lp - v, w, wp, p);
      }

      template<size_t N, enum ImplSelector Impl>
      void ntt(uint64_t *x, const uint64_t *w, const uint64_t *wp, const uint64_t p) {
        const uint64_t Lp = 2 * p;
        size_t t = N / 2;
        ++w, ++wp;

        for(size_t m = 1; m < N / 4; m *= 2) {
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + t;
          for(size_t i = 0; i < m; ++i) {
            for(size_t j = 0; j < t; j += 4) { // unrolling
              forward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              forward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              forward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              forward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
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
            forward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
            forward_butterfly_lazy<Impl>(xp0, xp1, *w, *wp, Lp, p);

            xp0 += 2 + 1, xp1 += 2 + 1; // combine the last xp0++, xp1++ to this increment
            ++w, ++wp;
          }
        }

        { // last loop, t = 1
          const size_t m = N / 2;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + 1;
          for(size_t i = 0; i < m; ++i) {
            forward_butterfly_lazy<Impl>(xp0, xp1, *w, *wp, Lp, p);

            xp0 += 2, xp1 += 2; // combine the last xp0++, xp1++ to this increment
            ++w, ++wp;
          }
        }
      }

      template<size_t N, enum ImplSelector Impl>
      void intt(uint64_t *x, const uint64_t *w, const uint64_t *wp,
                const uint64_t ni, const uint64_t nip, const uint64_t p) {

        const uint64_t Lp = 2 * p;

        ++w, ++wp;

        { // first loop, t = 1
          constexpr size_t m = N / 2;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + 1;
          for(size_t i = 0; i < m; ++i) {
            backward_butterfly_lazy<Impl>(xp0, xp1, *w, *wp, Lp, p);

            xp0 += 2, xp1 += 2;
            ++w, ++wp;
          }
        }

        { // second loop, t = 2
          constexpr size_t m = N / 4;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + 2;
          for(size_t i = 0; i < m; ++i) {
            backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
            backward_butterfly_lazy<Impl>(xp0, xp1, *w, *wp, Lp, p);

            xp0 += 2 + 1, xp1 += 2 + 1;
            ++w, ++wp;
          }
        }

        { // main loop, 4 <= t < n / 2
          size_t t = 4;
          for(size_t m = N / 8; m > 1; m /= 2) {
            uint64_t *xp0 = x;
            uint64_t *xp1 = x + t;
            for(size_t i = 0; i < m; ++i) {
              for(size_t j = 0; j < t; j += 4) {
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2;
          }
        }

        { // last loop t = n / 2,
          constexpr size_t t = N / 2;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + t;
          for(size_t j = 0; j < t; j += 4) {
            backward_butterfly_last<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
            backward_butterfly_last<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
            backward_butterfly_last<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
            backward_butterfly_last<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
          }
        }

      }

    }
  }
}

#endif //FASTERNTT_NTT_LN16_HPP
