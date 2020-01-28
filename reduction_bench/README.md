# Benchmark for Reduction Algorithms

This is a benchmark test for 64-bit x 64-bit modular multiplications. We consider following algorithms:
- Pseudo-Mersenne reduction[MOV96]
- Montgomery multiplication[MOV96]
- Barrett reduction[MOV96]
- NFLlib reduction[MBGG+16]
- K-REDC[LN16]
- Shoup's trick[Har14] (for fixed argument multiplication)
- gcc native (just for reference)

[MOV96] A. Menezes, P. van Oorschot, S. Vanstone, "Handbook of Applied Cryptography," CRC Press, 1996.  
[Har14] D. Harvey, "Faster arithmetic for number-theoretic transforms," Journal of Symbolic Computations, Vol.60, pp.113-119, 2014.  
[LN16] P. Longa and M. Naehrig, "Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography," CANS 2016, LNCS 10052, pp.124-139, 2016.  
[MBGG+16] C. A. Melchor, J. Barrier, S. Guelton, A. Guinet, M. Killijian, and T. Lepoint, "NFLlib: NTT-Based Fast Lattice Library," CT-RSA 2016, LNCS 9610, pp.341-356, 2016.

## Prerequisites

- [GMP](https://gmplib.org/)

## Results

All results include a 64-bit x 64-bit multiplication (which consumes almost 3 clks), 
and exclude conditional subtractions for outputs.

|Alg.|clk|
|----|----|
|Pseudo-Mersenne|18.733|
|Montgomery|**11.450**|
|Barrett|20.628|
|NFLlib|13.203|
|K-REDC|**7.539**|
|Shoup|**5.043**|
|gcc|92.432|

## Contact
Takuya HAYASHI (t-hayashi@nict.go.jp)
