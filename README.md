# FasterNTT - A Fast Implementation of Number-theoretic Transform

This is a fast implementation of NTT for Ring-LWE homomorphic encryption, FasterNTT, 
which is almost twice faster than NFLlib. FasterNTT employs

- Shoup's trick (for NTT), Montgomery multiplication (for point-wise mult.),
- and Scott's slothful NTT instead of Longa-Naehrig's NTT, which is employed in most of all (open-sourced) implementations.  

## Prerequisites

- [OpenSSL](https://www.openssl.org/)

## How to Use

This is a header only library, so just include and enjoy!

For benchmarks and tests, install following libraries:
- [GMP](https://gmplib.org/) (for test)
- [NTL](https://www.shoup.net/ntl/) (for test)
- [google/benchmark](https://github.com/google/benchmark) (for benchmark)
- [NFLlib](https://github.com/quarkslab/NFLlib) (if you want to compare)
- [YELL](https://github.com/fionser/YELL) (if you want to compare)

then build as usual cmake projects:
```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
```


## Abstract

Polynomial multiplication is the main bottleneck in Ring-LWE Cryptography, and most of all implementation researches focus on the Number-theoretic Transform (NTT) for efficient polynomial multiplications. The NTT consists of:

- modular arithmetic and
- butterfly computation.

Their performance has much influence on the performance of polynomial multiplication and so Ring-LWE (homomorphic) encryptions.

In this work, we comprehensively implemented efficient algorithms of modular arithmetic and NTTs in the Ring-LWE homomorphic encryptions case, and we found the combination of Shoup's trick (for NTT) and Montgomery multiplication (for point-wise multiplication) is the fastest way to compute Ring-LWE polynomial multiplication. Also, we found the Scott's slothful NTT is quite faster than the Longa-Naehrig's NTT, which most of all (open-sourced) implementations employ. Our experiment shows that our NTT implementation, which combines (modified) Scott's NTT, Shoup's trick, and Montgomery multiplication, is almost x2 faster than NFLlib, and x1.8 faster than YELL.

[[Har14]](https://doi.org/10.1016/j.jsc.2013.09.002) D. Harvey, "Faster arithmetic for number-theoretic transforms," Journal of Symbolic Computations, Vol.60, pp.113-119, 2014.  
[[LN16]](https://doi.org/10.1007/978-3-319-48965-0_8) P. Longa and M. Naehrig, "Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography," CANS 2016, LNCS 10052, pp.124-139, 2016.  
[[Sco17]](https://doi.org/10.1007/978-3-319-71045-7_13) M. Scott, "A Note on the Implementation of the Number Theoretic Transform," IMACC 2017, LNCS 10655, pp.247-258, 2017.

## Current result (n = 8192, #q = 4)

### Comparison with NFLlib and YELL

On Core i9-9980XE (3.00 GHz), turbo boost disabled, single thread, compiled by gcc 8.3.0:

|    | mean |
|----|----|
|FasterNTT NTT | **259 us** |
|FasterNTT INTT | **285 us** |
|FasterNTT PMUL | **36.1 us** |
|NFLlib NTT | 539 us |
|NFLlib INTT | 584 us |
|NFLlib PMUL | 55.5 us |
|YELL NTT | 474 us |
|YELL INTT | 468 us |
|YELL PMUL | 166 us |

On MacBook Pro 2017 (Core i7, 3.50 GHz), turbo boost enabled, single thread, compiled by gcc 9.2.0:

|    | mean |
|----|----|
|FasterNTT NTT | **204 us** |
|FasterNTT INTT | **226 us** |
|FasterNTT PMUL | **30.7 us** |
|NFLlib NTT |  354 us |
|NFLlib INTT | 404 us |
|NFLlib PMUL | 59.6 us |
|YELL NTT | 366 us |
|YELL INTT | 365 us |
|YELL PMUL | 125 us |

Note:

- NFLlib is compiled with `-DCMAKE_BUILD_TYPE=Release -DNFL_OPTIMIZED=ON`.
- YELL is compiled with `-DCMAKE_BUILD_TYPE=Release`.
- FasterNTT uses 59-bit q, whereas NFLlib and YELL use 62-bit q, which means FasterNTT naturally loses around 5% performance compared to NFLlib and YELL, but considering the fact ours is still about twice faster.

### Comparison with Several NTT Algorithms and Implementations

We implemented

- Longa-Naehrig's NTT (Longa-Naehrig's NTT + Shoup's trick + NFLlib reduction)
- Scott's slothful NTT (Scott's NTT + Shoup's trick + NFLlib reduction)
- FasterNTT ((modified) Scott's NTT + Shoup's trick + Montgomery multiplication)

using several optimized levels:

- Generic: using c++ (+ GCC extensions) only,
- x86_64: using inline assembler for x86_64 arch.,
- SSE: using SSE instructions,
- AVX2: using AVX2 instructions.

Note: SSE and AVX2 are implemented on FasterNTT only since it is the fastest. Also, point-wise multiplication using SIMD (SSE or AVX2) is not available because of their slowness (we implemented but significantly slower than x86_64 ver.).

The benchmark below is measured on Core i9-9980XE (3.00 GHz), turbo boost disabled, single thread, compiled by gcc 8.3.0.

| Alg. | Impl. |    | mean |
|----|----|----|----|
|Longa-Naehrig | Generic | NTT | 442 us |
|Longa-Naehrig | Generic | INTT | 558 us |
|Longa-Naehrig | Generic | PMUL | 55.4 us  |
|Longa-Naehrig | x86_64 | NTT | 345 us |
|Longa-Naehrig | x86_64 | INTT | 349 us |
|Longa-Naehrig | x86_64 | PMUL | 70.0 us  |
|Scott | Generic | NTT | 371 us |
|Scott | Generic | INTT | 378 us |
|Scott | Generic | PMUL | 55.3 us  |
|Scott | x86_64 | NTT | 314 us |
|Scott | x86_64 | INTT | 327 us |
|Scott | x86_64 | PMUL | 64.9 us  |
|FasterNTT | Generic | NTT | 371 us |
|FasterNTT | Generic | INTT | 362 us |
|FasterNTT | Generic | PMUL | **34.9 us** |
|FasterNTT | x86_64 | NTT | 317 us |
|FasterNTT | x86_64 | INTT | 313 us |
|FasterNTT | x86_64 | PMUL | **35.7 us** |
|FasterNTT | SSE | NTT | 448 us |
|FasterNTT | SSE | INTT | 486 us |
|FasterNTT | AVX2 | NTT | **260 us** |
|FasterNTT | AVX2 | INTT | **283 us** |

## Contact
Takuya HAYASHI (t-hayashi@nict.go.jp)
