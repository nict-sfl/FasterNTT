/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_ARITH_ASM_HPP
#define FASTERNTT_ARITH_ASM_HPP

namespace fntt {
  namespace core {
    // return (a >= b) ? a - b : a
    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t sub_conditional(uint64_t a, uint64_t b) {
      uint64_t t;
      __asm__ volatile (
      "mov %[a], %[t]\n\t"
      "sub %[b], %[t]\n\t"
      "cmovnc %[t], %[a]\n\t"
      : [a] "+&r"(a), [t] "=&r"(t)
      : [b] "r"(b)
      : "cc");
      return a;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t mulhi(uint64_t a, uint64_t b) {
      uint64_t t;
      __asm__ volatile (
        "mulx %[b], %[t], %[a]\n\t"
      : [a] "+&d"(a), [t] "=&r"(t)
      : [b] "r"(b)
      : );
      return a;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t redc_barrett(uint64_t a, const uint64_t r, const uint64_t p) {
      uint64_t t0, t1;
      __asm__ volatile (
      "mulx %[r], %[t0], %[t1]\n\t"
      "imul %[p], %[t1]\n\t"
      "sub %[t1], %[a]\n\t"
      "mov %[a], %[t0]\n\t"
      "sub %[p], %[t0]\n\t"
      "cmovnc %[t0], %[a]\n\t"
      : [a] "+&d"(a), [t0] "=&r"(t0), [t1] "=&r"(t1)
      : [r] "r"(r), [p] "r"(p)
      : "cc");
      return a;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t redc_barrett_lazy(uint64_t a, const uint64_t r, const uint64_t p) {
      uint64_t t0, t1;
      __asm__ volatile (
      "mulx %[r], %[t0], %[t1]\n\t"
      "imul %[p], %[t1]\n\t"
      "sub %[t1], %[a]\n\t"
      : [a] "+&d"(a), [t0] "=&r"(t0), [t1] "=&r"(t1)
      : [r] "r"(r), [p] "r"(p)
      : "cc");
      return a;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t mul_barrett_fixed(uint64_t a, const uint64_t b, const uint64_t bp,
                                      const uint64_t p) {
      uint64_t t0, t1;
      __asm__ volatile (
      "mulx %[bp], %[t0], %[t1]\n\t"
      "imul %[b], %[a]\n\t"
      "imul %[p], %[t1]\n\t"
      "sub %[t1], %[a]\n\t"
      "mov %[a], %[t0]\n\t"
      "sub %[p], %[t0]\n\t"
      "cmovnc %[t0], %[a]\n\t"
      : [a] "+&d"(a), [t0] "=&r"(t0), [t1] "=&r"(t1)
      : [b] "r"(b), [bp] "r"(bp), [p] "r"(p)
      : "cc");
      return a;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t mul_barrett_fixed_lazy(uint64_t a, const uint64_t b, const uint64_t bp, const uint64_t p) {
      uint64_t t0, t1;
      __asm__ volatile (
      "mulx %[bp], %[t0], %[t1]\n\t"
      "imul %[b], %[a]\n\t"
      "imul %[p], %[t1]\n\t"
      "sub %[t1], %[a]\n\t"
      : [a] "+&d"(a), [t0] "=&r"(t0), [t1] "=&r"(t1)
      : [b] "r"(b), [bp] "r"(bp), [p] "r"(p)
      : "cc");
      return a;
    }

    template<enum ImplSelector Impl, std::enable_if_t<Impl == x86_64 || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t mul_montgomery(const uint64_t a, const uint64_t b, const uint64_t m,
                                   const uint64_t p) {
      uint64_t t0, t1, c0, c1;
      __asm__ volatile (
      "mulx %[b], %[t0], %[t1]\n\t"
      "mov %[m], %[c0]\n\t"
      "mov %[p], %%rdx\n\t"
      "imul %[t0], %[c0]\n\t"
      "mulx %[c0], %[c0], %[c1]\n\t"
      "add %[t0], %[c0]\n\t"
      "adc %[t1], %[c1]\n\t"
      "mov %[c1], %[t0]\n\t"
      "sub %[p], %[t0]\n\t"
      "cmovnc %[t0], %[c1]\n\t"
      : [c1] "=&r"(c1), [c0] "=&r"(c0), [t1] "=&r"(t1), [t0] "=&r"(t0)
      : [a] "d"(a), [b] "r"(b), [m] "r"(m), [p] "r"(p)
      : "cc");
      return c1;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t mul_montgomery_lazy(const uint64_t a, const uint64_t b, const uint64_t m, const uint64_t p) {
      uint64_t t0, t1, c0, c1;
      __asm__ volatile (
      "mulx %[b], %[t0], %[t1]\n\t"
      "mov %[m], %[c0]\n\t"
      "mov %[p], %%rdx\n\t"
      "imul %[t0], %[c0]\n\t"
      "mulx %[c0], %[c0], %[c1]\n\t"
      "add %[t0], %[c0]\n\t"
      "adc %[t1], %[c1]\n\t"
      : [c1] "=&r"(c1), [c0] "=&r"(c0), [t1] "=&r"(t1), [t0] "=&r"(t0)
      : [a] "d"(a), [b] "r"(b), [m] "r"(m), [p] "r"(p)
      : "cc");
      return c1;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t mul_nfl(const uint64_t a, const uint64_t b, const uint64_t v, const uint8_t s, const uint64_t p) {
      uint64_t c0, c1, t0, t1, r;
      __asm__ volatile (
      "mulx %[b], %[c0], %%rdx\n\t"
      "mov %[c0], %[r]\n\t"
      "mulx %[v], %[t0], %[t1]\n\t"
      "shld %[s], %[c0], %%rdx\n\t"
      "shl %[s], %[c0]\n\t"
      "add %[c0], %[t0]\n\t"
      "adc %%rdx, %[t1]\n\t"
      "imul %[p], %[t1]\n\t"
      "sub %[t1], %[r]\n\t"
      "mov %[r], %[c0]\n\t"
      "sub %[p], %[c0]\n\t"
      "cmovnc %[c0], %[r]\n\t"
      : [r] "=&r"(r), [t0] "=&r"(t0), [t1] "=&r"(t1), [c0] "=&r"(c0), [c1] "=&r"(c1)
      : [a] "d"(a), [b] "r"(b), [s]"c"(s), [v]"r"(v), [p]"r"(p)
      : "cc");
      return r;
    }

    template<enum ImplSelector Impl,
      std::enable_if_t<Impl == x86_64 || Impl == SSE || Impl == AVX2, std::nullptr_t> = nullptr>
    inline uint64_t mul_nfl_lazy(const uint64_t a, const uint64_t b, const uint64_t v, const uint8_t s, const uint64_t p) {
      uint64_t c0, c1, t0, t1, r;
      __asm__ volatile (
      "mulx %[b], %[c0], %%rdx\n\t"
      "mov %[c0], %[r]\n\t"
      "mulx %[v], %[t0], %[t1]\n\t"
      "shld %[s], %[c0], %%rdx\n\t"
      "shl %[s], %[c0]\n\t"
      "add %[c0], %[t0]\n\t"
      "adc %%rdx, %[t1]\n\t"
      "imul %[p], %[t1]\n\t"
      "sub %[t1], %[r]\n\t"
      : [r] "=&r"(r), [t0] "=&r"(t0), [t1] "=&r"(t1), [c0] "=&r"(c0), [c1] "=&r"(c1)
      : [a] "d"(a), [b] "r"(b), [s]"c"(s), [v]"r"(v), [p]"r"(p)
      : "cc");
      return r;
    }
  } // namespace core
} // namespace fntt

#endif //FASTERNTT_ARITH_ASM_HPP
