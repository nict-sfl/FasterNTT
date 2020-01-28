/*
 * @file const.hpp 
 * @author Takuya HAYASHI <t-hayashi@nict.go.jp>
 */

#ifndef FASTERNTT_CONST_HPP
#define FASTERNTT_CONST_HPP


namespace fntt{
  enum NTTSelector{
    LN16, S17, FNTT
  };

  enum ImplSelector{
    Generic, x86_64, SSE, AVX2
  };

  static constexpr size_t kWordSize = 64UL;
} // namespace fntt

#endif //FASTERNTT_CONST_HPP
