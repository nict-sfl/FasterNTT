/*
 * @file param.hpp 
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 */

#ifndef FASTERNTT_RQ_PARAM_HPP
#define FASTERNTT_RQ_PARAM_HPP

#include <cmath>

#include "utils.hpp"
#include "core/const.hpp"

#include "core/prime59.hpp"
#include "core/prime59_nfl.hpp"
#include "core/prime62_nfl.hpp"

#include "core/arith.hpp"
#include "core/arith_asm.hpp"
#include "number_theory.hpp"

namespace fntt {

  namespace {
    template<size_t N, enum ImplSelector Impl>
    void compute_twiddle_factor_table(uint64_t *tbl, const uint64_t w, const uint64_t p) {
      auto len = (size_t) std::log2((double) N);
      uint64_t t = w, wp = core::mul_barrett_fixed_prep(w, p);
      tbl[0] = 1UL;
      for (size_t i = 1; i < N; ++i) {
        tbl[utils::bit_reverse(i, len)] = t;
        t = core::mul_barrett_fixed<Impl>(t, w, wp, p);
      }
    }

    constexpr size_t select_s(const enum NTTSelector NTT) {
      size_t ps = 0;
      switch (NTT) {
        case LN16:
          ps = core::prime62nfl::kPrimeBitSize;
          break;
        case S17:
          ps = core::prime59nfl::kPrimeBitSize;
          break;
        default:
          assert(true);
      }
      return kWordSize - ps;
    }
  }

  template<size_t N, size_t qNum, enum NTTSelector NTT, enum ImplSelector Impl>
  class RqParam {
  public:
    static constexpr size_t s = select_s(NTT);
    uint64_t p[qNum];
    uint64_t v0[qNum]; // floor(2^128 / p) mod 2^64
    uint64_t rdp[qNum]; // floor(2^64 / p)
    uint64_t *w[qNum]; // power series of 2n-th root of unity w
    uint64_t *wp[qNum];
    uint64_t *wi[qNum]; // w[]^-1 mod p
    uint64_t *wip[qNum];
    uint64_t ni[qNum]; // N^-1 mod p
    uint64_t nip[qNum];

    template<enum NTTSelector S = NTT, std::enable_if_t<S == LN16, std::nullptr_t> = nullptr>
    RqParam() {
      static_assert(utils::is_power_of_2(N), "N must be a power of two");
      static_assert(qNum < core::prime62nfl::kMaxPrimeNum, "qNum < kMaxPrimeNumLimb");
      for (size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        this->p[qIdx] = core::prime62nfl::p[qIdx];
        this->v0[qIdx] = core::prime62nfl::v[qIdx];
        this->w[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wp[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wi[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wip[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
      }
      this->init();
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == S17, std::nullptr_t> = nullptr>
    RqParam() {
      static_assert(utils::is_power_of_2(N), "N must be a power of two");
      static_assert(qNum < core::prime59nfl::kMaxPrimeNum, "qNum < kMaxPrimeNumLimb");
      for (size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        this->p[qIdx] = core::prime59nfl::p[qIdx];
        this->v0[qIdx] = core::prime59nfl::v[qIdx];
        this->w[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wp[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wi[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wip[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
      }
      this->init();
    }

    void init() {
      for (size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        // (N * 2^64)^-1 mod p
        this->ni[qIdx] = number_theory::inv(N, this->p[qIdx]);
        this->nip[qIdx] = core::mul_barrett_fixed_prep(this->ni[qIdx], this->p[qIdx]);

        // barrett parameter, 2^64 / p
        this->rdp[qIdx] = core::mul_barrett_fixed_prep(1, this->p[qIdx]);

        uint64_t w_ = 0;
        if (!number_theory::find_2nth_root_of_unity(w_, N, this->p[qIdx]))
          throw std::runtime_error("could not find 2n-th root of unity");

        compute_twiddle_factor_table<N, Impl>(this->w[qIdx], w_, this->p[qIdx]);
        compute_twiddle_factor_table<N, Impl>(this->wi[qIdx], number_theory::inv(w_, this->p[qIdx]), this->p[qIdx]);

        for (size_t i = 1; i < N; ++i)
          this->wp[qIdx][i] = core::mul_barrett_fixed_prep(this->w[qIdx][i], this->p[qIdx]);

        // merge N^-1 mod p to the w^i's
        this->wi[qIdx][1] = core::mul_barrett_fixed<Impl>(this->wi[qIdx][1], this->ni[qIdx], this->nip[qIdx],
                                                          this->p[qIdx]);
        for (size_t i = 1; i < N; ++i)
          this->wip[qIdx][i] = core::mul_barrett_fixed_prep(this->wi[qIdx][i], this->p[qIdx]);

        // reordering wi so that the access pattern at inverse NTT is sequential.
        uint64_t tmp[N], *ptr = tmp + 1;
        for (size_t i = N / 2; i > 0; i /= 2)
          for (size_t j = i; j < i * 2; ++j)
            *ptr++ = this->wi[qIdx][j];
        for (size_t i = 1; i < N; ++i) this->wi[qIdx][i] = tmp[i];

        ptr = tmp + 1;
        for (size_t i = N / 2; i > 0; i /= 2)
          for (size_t j = i; j < i * 2; ++j)
            *ptr++ = this->wip[qIdx][j];
        for (size_t i = 1; i < N; ++i) this->wip[qIdx][i] = tmp[i];
      }
    }

    ~RqParam() {
      for (size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        free(this->w[qIdx]);
        free(this->wp[qIdx]);
        free(this->wi[qIdx]);
        free(this->wip[qIdx]);
      }
    }
  };

  // FasterNTT needs additional parameters
  template<size_t N, size_t qNum, enum ImplSelector Impl>
  class RqParam<N, qNum, FNTT, Impl> {
  public:
    static constexpr size_t s = kWordSize - core::prime59::kPrimeBitSize;
    uint64_t p[qNum];
    uint64_t m[qNum];
    uint64_t rdp[qNum];
    uint64_t *w[qNum];
    uint64_t *wp[qNum];
    uint64_t *wi[qNum];
    uint64_t *wip[qNum];
    uint64_t ni[qNum];
    uint64_t nip[qNum];
    uint64_t r[qNum];
    uint64_t rp[qNum];

    RqParam() {
      static_assert(utils::is_power_of_2(N), "N must be a power of two");
      static_assert(qNum < core::prime59::kMaxPrimeNumLimb, "qNum < kMaxPrimeNumLimb");
      for (size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        this->p[qIdx] = core::prime59::p[qIdx];
        this->w[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wp[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wi[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
        this->wip[qIdx] = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N);
      }
      this->init();
    }

    void init() {
      for (size_t qIdx = 0; qIdx < qNum; ++qIdx) {

        // barrett parameter, 2^64 / p
        this->rdp[qIdx] = core::mul_barrett_fixed_prep(1, this->p[qIdx]);

        // montgomery parameter, -p^-1 mod 2^64
        this->m[qIdx] = core::compute_montgomery_param(this->p[qIdx]);

        // 2^64 mod p
        //this->r[qIdx] = ((unsigned __int128) 1UL << 64U) % this->p[qIdx];
        this->r[qIdx] = (uint64_t) (-this->p[qIdx]) % this->p[qIdx];
        this->rp[qIdx] = core::mul_barrett_fixed_prep(this->r[qIdx], this->p[qIdx]);

        // (N * 2^64)^-1 mod p, (2^64)^-1 mod p part is for mapping back from Montgomery domain
        this->ni[qIdx] = number_theory::inv(((unsigned __int128) N << 64U) % this->p[qIdx], this->p[qIdx]);
        this->nip[qIdx] = core::mul_barrett_fixed_prep(this->ni[qIdx], this->p[qIdx]);

        uint64_t w_ = 0;
        if (!number_theory::find_2nth_root_of_unity(w_, N, this->p[qIdx]))
          throw std::runtime_error("could not find 2n-th root of unity");

        compute_twiddle_factor_table<N, Impl>(this->w[qIdx], w_, this->p[qIdx]);
        compute_twiddle_factor_table<N, Impl>(this->wi[qIdx], number_theory::inv(w_, this->p[qIdx]), this->p[qIdx]);

        for (size_t i = 1; i < N / 2; ++i)
          this->wp[qIdx][i] = core::mul_barrett_fixed_prep(this->w[qIdx][i], this->p[qIdx]);
        // merge r mod p to the w's for mapping to Montgomery domain
        for (size_t i = N / 2; i < N; ++i) {
          this->w[qIdx][i] = core::mul_barrett_fixed<Impl>(this->w[qIdx][i], this->r[qIdx], this->rp[qIdx],
                                                           this->p[qIdx]);
          this->wp[qIdx][i] = core::mul_barrett_fixed_prep(this->w[qIdx][i], this->p[qIdx]);
        }

        // merge N^-1 mod p to the w^i's
        this->wi[qIdx][1] = core::mul_barrett_fixed<Impl>(this->wi[qIdx][1], this->ni[qIdx], this->nip[qIdx],
                                                          this->p[qIdx]);
        for (size_t i = 1; i < N; ++i)
          this->wip[qIdx][i] = core::mul_barrett_fixed_prep(this->wi[qIdx][i], this->p[qIdx]);

        // reordering wi so that the access pattern at inverse NTT is sequential.
        uint64_t tmp[N], *ptr = tmp + 1;
        for (size_t i = N / 2; i > 0; i /= 2)
          for (size_t j = i; j < i * 2; ++j)
            *ptr++ = this->wi[qIdx][j];
        for (size_t i = 1; i < N; ++i) this->wi[qIdx][i] = tmp[i];

        ptr = tmp + 1;
        for (size_t i = N / 2; i > 0; i /= 2)
          for (size_t j = i; j < i * 2; ++j)
            *ptr++ = this->wip[qIdx][j];
        for (size_t i = 1; i < N; ++i) this->wip[qIdx][i] = tmp[i];
      }
    }

    ~RqParam() {
      for (size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        free(this->w[qIdx]);
        free(this->wp[qIdx]);
        free(this->wi[qIdx]);
        free(this->wip[qIdx]);
      }
    }
  };

} // namespace fntt

#endif //FASTERNTT_RQ_PARAM_HPP
