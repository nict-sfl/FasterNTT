/*
 * @file ntt_core.hpp 
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 */

#ifndef FASTERNTT_NTT_HPP
#define FASTERNTT_NTT_HPP

#include "const.hpp"

#include "arith.hpp"
#include "arith_asm.hpp"

#include "arith_avx2.hpp"

namespace fntt {
  namespace core {
    namespace FNTT {

      template<enum ImplSelector Impl>
      inline void forward_butterfly_lazy(uint64_t *x, uint64_t *y,
                                         const uint64_t w, const uint64_t wp, const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = *x;
        v = mul_barrett_fixed_lazy<Impl>(*y, w, wp, p);
        *x = u + v;
        *y = u + Lp - v;
      }

      /// For the last layer:
      /// Twiddle factor w has been already multiplied by R, so u also needs to multiply by R,
      /// so that the output will be in Montgomery domain.
      template<enum ImplSelector Impl>
      inline void forward_butterfly_montgomery_domain_lazy(uint64_t *x, uint64_t *y,
                                                           uint64_t w, uint64_t wp, const uint64_t Lp,
                                                           const uint64_t r, const uint64_t rp, const uint64_t p) {
        uint64_t u, v;
        u = mul_barrett_fixed_lazy<Impl>(*x, r, rp, p);
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
        u = *x;
        v = *y;
        *x = redc_barrett_lazy<Impl>(u + v, r, p);
        *y = mul_barrett_fixed_lazy<Impl>(u + Lp - v, w, wp, p);
      }

      // Back from Montgomery domain, so (u + v) is multiplied by (NR)^-1 instead of N^-1.
      // Also the twiddle factor has been already multiplied by (NR)^-1.
      template<enum ImplSelector Impl>
      inline void backward_butterfly_montgomery_domain(uint64_t *x, uint64_t *y,
                                                       const uint64_t w, const uint64_t wp,
                                                       const uint64_t ni, const uint64_t nip,
                                                       const uint64_t Lp, const uint64_t p) {
        uint64_t u, v;
        u = *x;
        v = *y;
        *x = mul_barrett_fixed<Impl>(u + v, ni, nip, p);
        *y = mul_barrett_fixed<Impl>(u + Lp - v, w, wp, p);
      }

      /// Slothful NTT using Barrett reduction with fixed argument
      /// See: M. Scott, "A Note on Implementation of the Number Theoretic Transform", IMACC 2017.
      /// https://doi.org/10.1007/978-3-319-71045-7_13
      /// Output will be multiplied by R (means be in Montgomery domain) and in [0, 4q).
      template<size_t N, enum ImplSelector Impl, std::enable_if_t<
          Impl == Generic || Impl == x86_64, std::nullptr_t> = nullptr>
      void ntt(uint64_t *x, const uint64_t *w, const uint64_t *wp,
               const uint64_t r, const uint64_t rp, const uint64_t p) {

        const uint64_t Lp = 2 * p;
        size_t t = N / 2;

        ++w, ++wp;

        for (size_t m = 1; m < N / 4; m *= 2) {
          uint64_t *xp0 = x;
          uint64_t *xp1 = x + t;
          for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < t; j += 4) { // unrolling
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
          for (size_t i = 0; i < m; ++i) { // unrolling
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
          for (size_t i = 0; i < m; ++i) {
            forward_butterfly_montgomery_domain_lazy<Impl>(xp0, xp1, *w, *wp, Lp, r, rp, p);

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
      /// We slightly modified the algorithm so that the intermediate values are always storable in uint64_t.
      /// Input should be in Montgomery domain and in [0, Lq) where L = 2^(64 - 59 - 1) = 16
      /// Output will be fully reduced (be in [0, q)).
      template<size_t N, enum ImplSelector Impl, std::enable_if_t<
          Impl == Generic || Impl == x86_64, std::nullptr_t> = nullptr>
      void intt(uint64_t *x, const uint64_t *w, const uint64_t *wp,
                const uint64_t ni, const uint64_t nip, const uint64_t r, const uint64_t p) {

        constexpr size_t s = kWordSize - prime59::kPrimeBitSize - 1;
        constexpr size_t L = 1UL << s;
        const uint64_t Lp = L * p;
        uint64_t *xp0, *xp1;

        ++w, ++wp;

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
              for (size_t j = 0; j < t; j += 4) {
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          xp0 = x, xp1 = x + t;
          for (size_t i = 0; i < m; ++i) {
            size_t j = 0;
            for (; j < 2 * t / L; j += 4) {
              backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
              backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
              backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
              backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
            }
            for (; j < t; j += 4) {
              backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
            }
            xp0 += t, xp1 += t;
            ++w, ++wp;
          }
          t *= 2, m /= 2;

          for (int qq = 1; qq < Q; ++qq) {
            for (int ll = 0; ll < logL - 1; ++ll) {
              xp0 = x, xp1 = x + t;
              for (size_t i = 0; i < m; ++i) {
                size_t j = 0;
                for (; j < 1 * t / L; j += 4) {
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                }
                for (; j < 2 * t / L; j += 4) {
                    backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                    backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                    backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                    backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                }
                for (; j < t; j += 4) {
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                    backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                }
                xp0 += t, xp1 += t;
                ++w, ++wp;
              }
              t *= 2, m /= 2;
            }

            xp0 = x, xp1 = x + t;
            for (size_t i = 0; i < m; ++i) {
              size_t j = 0;
              for (; j < 2 * t / L; j += 4) {
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
              }
              for (; j < t; j += 4) {
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              }
              xp0 += t, xp1 += t;
              ++w, ++wp;
            }
            t *= 2, m /= 2;
          }

          for (int rr = 0; rr < R; ++rr) {
            xp0 = x, xp1 = x + t;
            for (size_t i = 0; i < m; ++i) {
              size_t j = 0;
              for (; j < 1 * t / L; j += 4) {
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
              }
              for (; j < 2 * t / L; j += 4) {
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
                  backward_butterfly_corr_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, r, p);
              }
              for (; j < t; j += 4) {
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
                  backward_butterfly_lazy<Impl>(xp0++, xp1++, *w, *wp, Lp, p);
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
            for (size_t j = 0; j < t; j += 4) { // unrolling
                backward_butterfly_montgomery_domain<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
                backward_butterfly_montgomery_domain<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
                backward_butterfly_montgomery_domain<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
                backward_butterfly_montgomery_domain<Impl>(xp0++, xp1++, *w, *wp, ni, nip, Lp, p);
            }
          }
        }
      }
    } // namespace FNTT
  } // namespace core
} // namespace fntt

#endif //FASTERNTT_NTT_HPP
