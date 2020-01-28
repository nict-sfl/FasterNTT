/*
 * @file uint64_arithmetic_asm.hpp 
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 */

#ifndef FASTERNTT_ARITH_HPP
#define FASTERNTT_ARITH_HPP

#include <cstdint>

namespace fntt {
  namespace core {

    // compute -(p^-1) mod 2^64
    uint64_t compute_montgomery_param(const uint64_t p) {
      uint64_t e = -(1UL << 63U) - 1, // Z/(2^64)Z^* has a multiplicative order of (2^64 - 2^63)
          a = -p, r = 1;
      while(e > 0) {
        if(e & 1U) r *= a;
        a *= a;
        e >>= 1U;
      }
      return r;
    }

    // compute floor(b * 2^64 / p)
    inline uint64_t mul_barrett_fixed_prep(const uint64_t b, const uint64_t p) {
      return ((unsigned __int128) b << 64U) / p;
    }

    typedef union {
      unsigned __int128 d128;
      uint64_t d64[2];
    } uint128_t;

    inline uint128_t mul128(const uint64_t a, const uint64_t b) {
      uint128_t r;
      r.d128 = (unsigned __int128) a * b;
      return r;
    }

    inline uint128_t adds128(uint128_t a, const uint128_t b) {
      a.d128 += b.d128;
      return a;
    }

    inline uint128_t lshift128(uint128_t a, const size_t s) {
      a.d128 <<= s;
      return a;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    inline uint64_t sub_conditional(const uint64_t a, const uint64_t b) {
      return (a >= b) ? a - b : a;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    inline uint64_t redc_barrett(const uint64_t a, const uint64_t r, const uint64_t p) {
      uint128_t t;
      uint64_t ret;
      t = mul128(a, r);
      ret = a - t.d64[1] * p;
      return (ret >= p) ? ret - p : ret;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    inline uint64_t redc_barrett_lazy(const uint64_t a, const uint64_t r, const uint64_t p) {
      uint128_t t;
      t = mul128(a, r);
      return a - t.d64[1] * p;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    inline uint64_t mul_barrett_fixed(const uint64_t a, const uint64_t b, const uint64_t bp, const uint64_t p) {
      uint128_t t;
      uint64_t r;
      t = mul128(a, bp);
      r = a * b - t.d64[1] * p;
      return (r >= p) ? r - p : r;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    inline uint64_t mul_barrett_fixed_lazy(const uint64_t a, const uint64_t b, const uint64_t bp, const uint64_t p) {
      uint128_t t;
      t = mul128(a, bp);
      return a * b - t.d64[1] * p;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    uint64_t mul_montgomery(const uint64_t a, const uint64_t b, const uint64_t m, const uint64_t p) {
      uint128_t t, c;
      uint64_t u, r;
      t = mul128(a, b);
      u = t.d64[0] * m;
      c = mul128(u, p);
      c = adds128(c, t);
      r = c.d64[1];
      return (r >= p) ? (r - p) : r;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    uint64_t mul_montgomery_lazy(const uint64_t a, const uint64_t b, const uint64_t m, const uint64_t p) {
      uint128_t t, c;
      uint64_t u;
      t = mul128(a, b);
      u = t.d64[0] * m;
      c = mul128(u, p);
      c = adds128(c, t);
      return c.d64[1];
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    uint64_t mul_nfl(const uint64_t a, const uint64_t b, const uint64_t v, const size_t s, const uint64_t p) {
      uint128_t c, t;
      uint64_t t0, t1, r;
      c = mul128(a, b);
      t0 = c.d64[0];
      t1 = c.d64[1];
      c = lshift128(c, s);
      t = mul128(v, t1);
      t = adds128(t, c);
      r = t0 - t.d64[1] * p;
      return (r >= p) ? r - p : r;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == Generic, std::nullptr_t> = nullptr>
    uint64_t mul_nfl_lazy(const uint64_t a, const uint64_t b, const uint64_t v, const size_t s, const uint64_t p) {
      uint128_t c, t;
      uint64_t t0, t1;
      c = mul128(a, b);
      t0 = c.d64[0];
      t1 = c.d64[1];
      c = lshift128(c, s);
      t = mul128(v, t1);
      t = adds128(t, c);
      return t0 - t.d64[1] * p;
    }
  } // namespace core
} // namespace fntt

#endif //FASTERNTT_ARITH_HPP
