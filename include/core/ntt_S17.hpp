/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_NTT_S17_HPP
#define FASTERNTT_NTT_S17_HPP

#include "const.hpp"

#ifdef USE_ASM
#include "arith_asm.hpp"
#else

#include "arith.hpp"

#endif

namespace fntt {
  namespace core {
    namespace S17 {

      template<enum ImplSelector Impl>
      inline void forward_butterfly_lazy(uint64_t *x, uint64_t *y,
                                         const uint64_t w, const uint64_t wp, const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = *x;
        v = mul_barrett_fixed_lazy<Impl>(*y, w, wp, p);
        *x = u + v;
        *y = u + Lp - v;
      }

      template<enum ImplSelector Impl>
      inline void forward_butterfly_lazy_last(uint64_t *x, uint64_t *y,
                                              const uint64_t w, const uint64_t wp, const uint64_t Lp,
                                              const uint64_t r, const uint64_t p) {
        uint64_t u, v;
        u = redc_barrett_lazy<Impl>(*x, r, p);
        v = mul_barrett_fixed_lazy<Impl>(*y, w, wp, p);
        *x = u + v;
        *y = u + Lp - v;
      }

      template<enum ImplSelector Impl>
      inline void backward_butterfly_lazy(uint64_t *x, uint64_t *y,
                                          const uint64_t w, const uint64_t wp,
                                          const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = *x;
        v = *y;
        *x = u + v;
        *y = mul_barrett_fixed_lazy<Impl>(u + Lp - v, w, wp, p);
      }

      // the u and v are reduced to [0, 2q) then compute butterfly.
      template<enum ImplSelector Impl>
      inline void backward_butterfly_corr_lazy(uint64_t *x, uint64_t *y,
                                               const uint64_t w, const uint64_t wp,
                                               const uint64_t Lp, const uint64_t r,
                                               const uint64_t p) {
        uint64_t u, v;
        u = redc_barrett_lazy<Impl>(*x, r, p);
        v = redc_barrett_lazy<Impl>(*y, r, p);
        *x = u + v;
        *y = mul_barrett_fixed_lazy<Impl>(u + Lp - v, w, wp, p);
      }

      template<enum ImplSelector Impl>
      inline void backward_butterfly(uint64_t *x, uint64_t *y,
                                     const uint64_t w, const uint64_t wp,
                                     const uint64_t ni, const uint64_t nip,
                                     const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = *x;
        v = *y;
        *x = mul_barrett_fixed<Impl>(u + v, ni, nip, p);
        *y = mul_barrett_fixed<Impl>(u + Lp - v, w, wp, p);
      }

      // for the last layer
      template<enum ImplSelector Impl>
      inline void backward_butterfly_corr(uint64_t *x, uint64_t *y,
                                          const uint64_t w, const uint64_t wp,
                                          const uint64_t ni, const uint64_t nip,
                                          const uint64_t Lp, const uint64_t r,
                                          const uint64_t p) {
        uint64_t u, v;
        u = redc_barrett_lazy<Impl>(*x, r, p);
        v = redc_barrett_lazy<Impl>(*y, r, p);
        *x = mul_barrett_fixed<Impl>(u + v, ni, nip, p);
        *y = mul_barrett_fixed<Impl>(u + Lp - v, w, wp, p);
      }

      /// Slothful NTT using Barrett reduction with fixed argument
      /// See: M. Scott, "A Note on Implementation of the Number Theoretic Transform", IMACC 2017.
      /// https://doi.org/10.1007/978-3-319-71045-7_13
      /// Output will be in [0, 4p).
      /// Suppose (2 * (log2(N/2) + 1) * p <= 2^64 so that carries during NTT can be stored in uint64_t w/o reductions.
      template<size_t N, enum ImplSelector Impl>
      void ntt(uint64_t *x, const uint64_t *w, const uint64_t *wp, const uint64_t rdp, const uint64_t p) {

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

        { // last loop, t = 1, each element is in [0, (2 * log2(N / 2) + 1) * p)
          const size_t m = N / 2;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + 1;
          for(size_t i = 0; i < m; ++i) {
            forward_butterfly_lazy_last<Impl>(xp0, xp1, *w, *wp, Lp, rdp, p);

            xp0 += 2, xp1 += 2; // combine the last xp0++, xp1++ to this increment
            ++w, ++wp;
          }
        }
      }

      /// (Semi-)slothful inverse NTT using Barrett reduction with fixed argument
      /// Since inverse NTT needs more capacity to store intermediate values during lazy reductions,
      /// it is necessary to compute additional reductions, and so we add the prefix "semi-".
      /// See: M. Scott, "A Note on Implementation of the Number Theoretic Transform", IMACC 2017.
      /// https://doi.org/10.1007/978-3-319-71045-7_13
      template<size_t N, enum ImplSelector Impl>
      void intt(uint64_t *x, const uint64_t *w, const uint64_t *wp, const uint64_t ni, const uint64_t nip,
                const uint64_t r, const uint64_t p) {

        constexpr size_t L = 1UL << (kWordSize - prime59nfl::kPrimeBitSize - 1);
        const uint64_t Lp = L * p;

        ++w, ++wp;

        { // first loop, t = 1
          constexpr size_t m = N / 2;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + 1;
          for(size_t i = 0; i < m; ++i) {
            // backward_butterfly_lazy<Impl>(xp0, xp1, *w, *wp, Lp, p);
            // Input will be in [0, Lq) so must be reduced at first
            backward_butterfly_corr_lazy<Impl>(xp0, xp1, *w, *wp, Lp, r, p);
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
          size_t m = N / 8;

          for(; m >= (N / L); m /= 2) { // there's no need of additional reductions
            uint64_t *xp0 = x;
            uint64_t *xp1 = x + t;
            for(size_t i = 0; i < m; ++i) {
              for(size_t j = 0; j < t; j += 4) { // unrolling
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

          for(; m > 1; m /= 2) { // some needs reductions, some don't.
            uint64_t *xp0 = x;
            uint64_t *xp1 = x + t;
            for(size_t i = 0; i < m; ++i) {
              size_t j = 0;
              for(; j < ((N / L) / (2 * m)); j += 4) { // unrolling
                backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
              }
              for(; j < t; j += 4) { // unrolling
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

        { // last loop t = n / 2, some needs reductions, some don't.
          constexpr size_t m = 1;
          constexpr size_t t = N / 2;
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + t;

          size_t j = 0;
          for(; j < ((N / L) / (2 * m)); j += 4) { // unrolling
            backward_butterfly_corr<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, r, p);
            backward_butterfly_corr<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, r, p);
            backward_butterfly_corr<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, r, p);
            backward_butterfly_corr<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, r, p);
          }
          for(; j < t; j += 4) { // unrolling
            backward_butterfly<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
            backward_butterfly<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
            backward_butterfly<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
            backward_butterfly<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
          }
        }
      }

    }
  } // namespace core
} // namespace fntt

#endif //FASTERNTT_NTT_S17_HPP
