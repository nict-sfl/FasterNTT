/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_UTILS_SIMD_HPP
#define FASTERNTT_UTILS_SIMD_HPP

#include <immintrin.h>

inline void print(const __m256i a) {
  uint64_t t[4];
  _mm256_storeu_si256((__m256i *) t, a);
  printf("%lx %lx %lx %lx\n", t[0], t[1], t[2], t[3]);
}

inline int shuffle_imm(const unsigned i0, const unsigned i1, const unsigned i2, const unsigned i3){
  return (int)((i0 << 0U) | (i1 << 2U) | (i2 << 4U) | (i3 << 6U));
}

inline int blend32_imm(const unsigned i0, const unsigned i1, const unsigned i2, const unsigned i3,
                       const unsigned i4, const unsigned i5, const unsigned i6, const unsigned i7) {
  return (int)((i0 << 0U) | (i1 << 1U) | (i2 << 2U) | (i3 << 3U) | (i4 << 4U) | (i5 << 5U) | (i6 << 6U) | (i7 << 7U));
}

#endif //FASTERNTT_UTILS_SIMD_HPP
