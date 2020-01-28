#include "test.h"

int main() {
  const size_t niter = 1 << 24;

  test_pm(niter);
  test_pm2(niter);
  test_mm(niter);
  test_brt(niter);
  test_brt2(niter);
  test_brt_fixed(niter);
  test_nfl(niter);
  test_kred(niter);

  test_brt_fixed_asm(niter);
  test_mm_asm(niter);
  test_nfl_asm(niter);
  test_kred_asm(niter);

  return 0;
}