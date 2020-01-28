/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_RQ_HPP
#define FASTERNTT_RQ_HPP

#include <cstdint>
#include <cmath>
#include <cassert>
#include <cstring>
#include <stdexcept>
#include <random>

#include "utils.hpp"
#include "core/const.hpp"

#include "csprng.hpp"
#include "Rq_param.hpp"

#include "core/arith.hpp"
#include "core/arith_asm.hpp"

#include "core/ntt_LN16.hpp"
#include "core/ntt_S17.hpp"
#include "core/fntt.hpp"
#include "core/fntt_sse.hpp"
#include "core/fntt_avx2.hpp"
#include "core/pmul.hpp"

namespace fntt {
  template<size_t N, size_t qNum, enum NTTSelector NTT = FNTT, enum ImplSelector Impl = AVX2>
  class Rq {
  public:
    static csprng::sampler<csprng::aes128> sampler_;
    static RqParam<N, qNum, NTT, Impl> param_;
    uint64_t *data_;

    Rq() : data_(nullptr) {
      this->data_ = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N * qNum);
    }

    Rq(const Rq &a) : data_(nullptr) {
      this->data_ = (uint64_t *) utils::aligned_malloc(utils::AVX_MEMORY_ALIGNMENT, sizeof(uint64_t) * N * qNum);
      memcpy(this->data_, a.data_, sizeof(uint64_t) * N * qNum);
    }

    ~Rq() {
      free(this->data_);
    }

    void cleanse() {
      memset(this->data_, sizeof(uint64_t) * N * qNum, 0);
      // memset_s(this->data_, sizeof(uint64_t) * N * qNum, 0, sizeof(uint64_t) * N * qNum);
    }

    void set(const Rq &b) {
      memcpy(this->data_, b.data_, sizeof(uint64_t) * N * qNum);
    }

    void print(const size_t num = N) const {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        for(size_t i = 0; i < num; ++i) printf("%lu %5lu: %lu\n", qIdx, i, (*this)(qIdx, i));
      }
    }

    Rq &operator=(const Rq &b) {
      if(this != &b) this->set(b);
      return *this;
    }

    bool operator==(const Rq &b) const {
      return cmp(*this, b);
    }

    bool operator!=(const Rq &b) const {
      return !cmp(*this, b);
    }

    uint64_t *operator[](const size_t i) {
      return this->data_ + i * N;
    }

    const uint64_t *operator[](const size_t i) const {
      return this->data_ + i * N;
    }

    uint64_t &operator()(const size_t i, const size_t j) {
      return this->operator[](i)[j];
    }

    uint64_t operator()(const size_t i, const size_t j) const {
      return this->operator[](i)[j];
    }

    static bool cmp(const Rq &a, const Rq &b) {
      bool r = true;
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        for(size_t i = 0; i < N; ++i) r &= a(qIdx, i) == b(qIdx, i);
      }
      return r;
    }

    constexpr unsigned int deg() const {
      return N;
    }

    void rand() {
      this->sampler_.operator()((unsigned char*)this->data_, sizeof(uint64_t) * N * qNum);
      this->redc();
    }

    static void init() {
      Rq::param_.init();
    }

    void redc() {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        uint64_t *ptr = (*this)[qIdx];
        for(size_t i = 0; i < N; ++i)
          this->data_[i + N * qIdx] = core::redc_barrett<Impl>(this->data_[i + N * qIdx],
                                                               this->param_.rdp[qIdx], this->param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == LN16, std::nullptr_t> = nullptr>
    void ntt() {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::LN16::ntt<N, Impl>(this->data_ + N * qIdx, this->param_.w[qIdx], this->param_.wp[qIdx],
                                 this->param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == S17, std::nullptr_t> = nullptr>
    void ntt() {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::S17::ntt<N, Impl>(this->data_ + N * qIdx, this->param_.w[qIdx], this->param_.wp[qIdx],
                                this->param_.rdp[qIdx], this->param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == FNTT, std::nullptr_t> = nullptr>
    void ntt() {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::FNTT::ntt<N, Impl>(this->data_ + N * qIdx, this->param_.w[qIdx], this->param_.wp[qIdx],
                                 this->param_.r[qIdx], this->param_.rp[qIdx], this->param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == LN16, std::nullptr_t> = nullptr>
    void intt() {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::LN16::intt<N, Impl>(this->data_ + N * qIdx, this->param_.wi[qIdx], this->param_.wip[qIdx],
                                  this->param_.ni[qIdx], this->param_.nip[qIdx], this->param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == S17, std::nullptr_t> = nullptr>
    void intt() {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::S17::intt<N, Impl>(this->data_ + N * qIdx, this->param_.wi[qIdx], this->param_.wip[qIdx],
                                 this->param_.ni[qIdx], this->param_.nip[qIdx], this->param_.rdp[qIdx],
                                 this->param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == FNTT, std::nullptr_t> = nullptr>
    void intt() {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::FNTT::intt<N, Impl>(this->data_ + N * qIdx, this->param_.wi[qIdx], this->param_.wip[qIdx],
                                  this->param_.ni[qIdx], this->param_.nip[qIdx], this->param_.rdp[qIdx],
                                  this->param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == LN16 || S == S17, std::nullptr_t> = nullptr>
    static void pmul(Rq &c, const Rq &a, const Rq &b) {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::pmul_nfl<N, Impl>(c.data_ + N * qIdx, a.data_ + N * qIdx, b.data_ + N * qIdx,
                                Rq::param_.v0[qIdx], Rq::param_.s, Rq::param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == FNTT, std::nullptr_t> = nullptr>
    static void pmul(Rq &c, const Rq &a, const Rq &b) {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::pmul_montgomery<N, Impl>(c.data_ + N * qIdx, a.data_ + N * qIdx, b.data_ + N * qIdx,
                                       Rq::param_.m[qIdx], Rq::param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == LN16 || S == S17, std::nullptr_t> = nullptr>
    static void pmul_lazy(Rq &c, const Rq &a, const Rq &b) {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::pmul_lazy_nfl<N, Impl>(c.data_ + N * qIdx, a.data_ + N * qIdx, b.data_ + N * qIdx,
                                     Rq::param_.v0[qIdx], Rq::param_.s, Rq::param_.p[qIdx]);
      }
    }

    template<enum NTTSelector S = NTT, std::enable_if_t<S == FNTT, std::nullptr_t> = nullptr>
    static void pmul_lazy(Rq &c, const Rq &a, const Rq &b) {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        core::pmul_montgomery_lazy<N, Impl>(c.data_ + N * qIdx, a.data_ + N * qIdx, b.data_ + N * qIdx,
                                            Rq::param_.m[qIdx], Rq::param_.p[qIdx]);
      }
    }

    static void mul(Rq &c, const Rq &a, const Rq &b) {
      Rq a_(a), b_(b);
      a_.ntt(), b_.ntt(); // in [0, 4q)
      pmul_lazy(c, a_, b_); // in [0, 2q)
      c.intt(); // in [0, q)
    }

    static void add(Rq &c, const Rq &a, const Rq &b) {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        for(size_t i = 0; i < N; ++i) {
          c.data_[i + N * qIdx] =
            core::redc_barrett<Impl>(a.data_[i + N * qIdx] + b.data_[i + N * qIdx], Rq::param_.rdp[qIdx],
                                     Rq::param_.p[qIdx]);
        }
      }
    }

    static void add_lazy(Rq &c, const Rq &a, const Rq &b) {
      for(size_t i = 0; i < N * qNum; ++i) c.data_[i] = a.data_[i] + b.data_[i];
    }

    // @todo: optimize sub, sub_lazy, determine Lp
    static void sub(Rq &c, const Rq &a, const Rq &b) {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        const uint64_t Lp = 2 * Rq::param_.p[qIdx];
        for(size_t i = 0; i < N; ++i) {
          const_cast<Rq &>(b).data_[i + N * qIdx] =
            core::redc_barrett_lazy<Impl>(b.data_[i + N * qIdx], Rq::param_.rdp[qIdx], Rq::param_.p[qIdx]);
          c.data_[i + N * qIdx] = core::redc_barrett<Impl>(a.data_[i + N * qIdx] + Lp - b.data_[i + N * qIdx],
                                                           Rq::param_.rdp[qIdx], Rq::param_.p[qIdx]);
        }
      }
    }

    static void sub_lazy(Rq &c, const Rq &a, const Rq &b) {
      for(size_t qIdx = 0; qIdx < qNum; ++qIdx) {
        const uint64_t Lp = 2 * Rq::param_.p[qIdx];
        for(size_t i = 0; i < N; ++i) {
          const_cast<Rq &>(b).data_[i + N * qIdx] =
            core::redc_barrett_lazy<Impl>(b.data_[i + N * qIdx], Rq::param_.rdp[qIdx], Rq::param_.p[qIdx]);
          c.data_[i + N * qIdx] = a.data_[i + N * qIdx] + Lp - b.data_[i + N * qIdx];
        }
      }
    }
  };

  template<size_t N, size_t qNum, enum NTTSelector NTT, enum ImplSelector Impl>
  csprng::sampler<csprng::aes128> Rq<N, qNum, NTT, Impl>::sampler_;
  template<size_t N, size_t qNum, enum NTTSelector NTT, enum ImplSelector Impl>
  RqParam<N, qNum, NTT, Impl> Rq<N, qNum, NTT, Impl>::param_;
} // namespace fntt

#endif //FASTERNTT_RQ_HPP
