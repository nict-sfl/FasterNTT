/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_UTILS_HPP
#define FASTERNTT_UTILS_HPP

namespace utils {
  static constexpr size_t AVX_MEMORY_ALIGNMENT = 32;

  void *aligned_malloc(const size_t alignment, const size_t size) {
#ifdef __USE_ISOC11
    void *ptr = aligned_alloc(alignment, size);
    if(ptr == nullptr) throw std::bad_alloc();
#else
    void *ptr = nullptr;
    if(posix_memalign(&ptr, alignment, size)) throw std::bad_alloc();
#endif
    return ptr;
  }

  // compute bit reverse of i of length len
  // ex.) bit_reverse(0x05, 8) = 0xa0
  size_t bit_reverse(const size_t i, const size_t len) {
    size_t v = i, r_ = i;
    size_t s = len - 1;

    for(v >>= 1U; v; v >>= 1U) {
      r_ <<= 1U;
      r_ |= v & 1U;
      s--;
    }
    return (r_ << s) & ((1UL << len) - 1);
  }

  constexpr bool is_power_of_2(size_t N){
    if(N == 0) return false;
    while(N > 1){
      if(N & 1U) return false;
      N >>= 1U;
    }
    return true;
  }

  constexpr size_t bitLength(size_t N) {
    if (N == 0) return 0;
    size_t l = 0;
    while (N > 1) {
      N >>= 1U;
      ++l;
    }
    return l;
  }

} // namespace utils

#endif //FASTERNTT_UTILS_HPP
