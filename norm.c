#include <stdio.h>
#include <math.h>

#include "norm.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

double inf_norm(double const * u, mesh const * m,
                double (*u_exact)(double, double))
{
  size_t i;
  double * coords;
  double max_err = 0, err;

  if (u == NULL || m == NULL)
    err_exit("No solution vector or mesh given!");

  for (i = 0; i < m->n_vertices; ++i)
  {
    coords = &m->coords[2*i];
    err = fabs(u_exact(coords[0], coords[1]) - u[i]);
    max_err = err > max_err ? err : max_err;
  }
  return max_err;
}

double l2_norm(double const * u, mesh const * m,
               double (*u_exact)(double, double))
{
  size_t i;
  double * coords;
  double err = 0., local_err;

  if (u == NULL || m == NULL)
    err_exit("No solution vector or mesh given!");

  for (i = 0; i < m->n_vertices; ++i)
  {
    coords = &m->coords[2*i];
    local_err = u_exact(coords[0], coords[1]) - u[i];
    err += local_err * local_err;
  }
  return sqrt(err)/(double)(m->n_vertices);
}
