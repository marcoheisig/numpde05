#include <stdio.h>
#include <math.h>

#include "exercise5.h"

/* Choose the testcase:
 *  1: u = x + y (zero RHS, should be solved up to machine precision)
 *  2: u = sin(pi*x)*sin(pi*y) (zero Dirichlet BC)
 *  3: u = cos(7*x)*cos(7*y) (non-zero DBC and RHS)
 */
#define TESTCASE 2

double u(double x, double y)
{
#if TESTCASE == 1
  return x + y;
#elif TESTCASE == 2
  return sin(M_PI * x) * sin(M_PI * y);
#elif TESTCASE == 3
  return cos(7. * x) * cos(7. * y);
#else
#error("Invalid value for TESTCASE!");
#endif
}

double f(double x, double y)
{
#if TESTCASE == 1
  return 0*(x + y);
#elif TESTCASE == 2
  return 2. * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
#elif TESTCASE == 3
  return 98. * cos(7. * x) * cos(7. * y);
#else
#error("Invalid value for TESTCASE!");
#endif
}

double g(unsigned char i, double x, double y)
{
#if TESTCASE == 1
  return x + y + 0*i;
#elif TESTCASE == 2
  return sin(M_PI * x) * sin(M_PI * y) + 0*i;
#elif TESTCASE == 3
  return cos(7. * x) * cos(7. * y) + 0*i;
#else
#error("Invalid value for TESTCASE!");
#endif
}

int main()
{
    size_t n_avail[] = { 3, 10, 20 , 40, 80 };
  size_t it;
  double errors[2];

  printf("+------+------------+------------+\n");
  printf("| n    | l2_norm    | inf_norm   |\n");
  printf("+------+------------+------------+\n");

  for (it = 0; it < sizeof(n_avail)/sizeof(size_t); ++it)
  {
    fem(n_avail[it], errors, &f, &g, &u);
    printf("| %4d | %8.4E | %8.4E |\n", (int)n_avail[it], errors[0], errors[1]);
  }

  printf("+------+------------+------------+\n");

  return 0;
}
