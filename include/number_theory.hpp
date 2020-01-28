/*
 * @file numth.hpp 
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 */

#ifndef FASTERNTT_NUMBER_THEORY_HPP
#define FASTERNTT_NUMBER_THEORY_HPP
namespace fntt {

  namespace number_theory {
    static constexpr uint64_t kMaxNiterRootFinding = 1UL << 20U;

    // compute a^e mod p by binary method
    uint64_t pow(uint64_t a, uint64_t e, const uint64_t p) {
      uint64_t r = 1;
      while(e > 0) {
        if(e & 1U) r = ((unsigned __int128) r * a) % p;
        a = ((unsigned __int128) a * a) % p;
        e >>= 1U;
      }
      return r;
    }
    // compute a^-1 mod p = a^(p - 2) mod p
    uint64_t inv(const uint64_t a, const uint64_t p) {
      return pow(a, p - 2, p);
    }

    bool is_2nth_root_of_unity(uint64_t a, const unsigned int n, const uint64_t p) {
      auto len = (size_t) std::log2((double) n);
      for(size_t i = 0; i <= len; ++i) {
        if(a == 1) return false;
        a = ((unsigned __int128) a * a) % p;
      }
      return a == 1;
    }

    bool find_2nth_root_of_unity(uint64_t &a, const unsigned int n, const uint64_t p) {
      const uint64_t e = (p - 1) / (2 * n);
      for(uint64_t a_ = 2; a_ < kMaxNiterRootFinding; ++a_) {
        a = a_;
        if(is_2nth_root_of_unity(a, n, p)) return true;
        a = pow(a, e, p);
        if(is_2nth_root_of_unity(a, n, p)) return true;
      }
      return false;
    }

  } // namespace number_theory

} // namespace fntt


#endif //FASTERNTT_NUMBER_THEORY_HPP
