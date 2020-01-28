#include "bench.h"

int main() {
  const size_t niter = 1 << 30;

  bench_mul128(niter);
  bench_gcc(niter);
  bench_kred(niter);
  bench_pm(niter);
  bench_pm2(niter);
  bench_nfl(niter);
  bench_mm(niter);
  bench_brt(niter);
  bench_brt2(niter);
  bench_brt_fixed(niter);

  bench_kred_asm(niter);
  bench_nfl_asm(niter);
  bench_mm_asm(niter);
  bench_mul_brt_fixed_asm(niter);

  return 0;
}