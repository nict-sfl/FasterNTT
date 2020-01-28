/**
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 **/

#ifndef FASTERNTT_CSPRNG_HPP
#define FASTERNTT_CSPRNG_HPP

// #define SAMPLER_SEED_STATIC

#include <random>
#include <cstring>
#include <openssl/evp.h>

namespace csprng {

  namespace {
    inline void set_seed(unsigned char *p, const size_t len) {
      std::random_device rd;
      for (size_t i = 0; i < len / 4; ++i) {
        *((uint32_t *) (p) + i) = (uint32_t) rd();
      }
#ifdef SAMPLER_SEED_STATIC
      memset(p, 0, len);
#endif
    }
  }

  class aes128{
    unsigned char internal_state[16] = {};
    EVP_CIPHER_CTX *ctx = nullptr;
  public:
    aes128(){
      ctx = EVP_CIPHER_CTX_new();
      set_seed(this->internal_state, 16);
      EVP_EncryptInit_ex(ctx, EVP_aes_128_ctr(), nullptr, this->internal_state, nullptr);
    }

    ~aes128(){
      EVP_CIPHER_CTX_free(ctx);
    }

    size_t size(){
      return 16;
    }

    size_t operator()(unsigned char *r){
      unsigned char zero[16] = {0};
      int olen;
      EVP_EncryptUpdate(ctx, r, &olen, zero, 16);
      return this->size();
    }
  };

  template <typename T>
  class sampler {
    T engine;
  public:
    void operator()(unsigned char *p, const size_t len) {
      const size_t size = this->engine.size();
      const size_t q = len / size, r = len % size;
      for(size_t i = 0; i < q; ++i) p += this->engine(p);
      if(r){
        unsigned char t[size];
        this->engine(t);
        memcpy(p, t, r);
      }
    }

  };
}

#endif //FASTERNTT_CSPRNG_HPP
