#include <stdio.h>
#include <string.h>

#include "fem.h"
#include "norm.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

void fem(size_t n, double errors[2], double (*fn_f)(double, double),
         double (*fn_g)(unsigned char, double, double),
         double (*fn_u)(double, double))
{
  mesh m;
  crs_matrix mat;
  double * u, * rhs;
  double local_stiffness[3][3];
  double local_load[3];
  size_t elem;

  /* 1. Allocate and generate mesh */
  get_mesh(&m, n);

#ifdef PRINT_DEBUG
  print_mesh(&m);
#endif

  /* 2. Allocate the linear system */
  init_matrix(&mat, &m);

  u = (double *) malloc(sizeof(double) * m.n_vertices);
  if (u == NULL) err_exit("Allocation of solution vector failed!");
  memset(u, 0, sizeof(double) * m.n_vertices);

  rhs = (double *) malloc(sizeof(double) * m.n_vertices);
  if (rhs == NULL)
    err_exit("Allocation of right hand side failed!");
  memset(rhs, 0, sizeof(double) * m.n_vertices);

  /* 3. Assemble the matrix */
  for (elem = 0; elem <  m.n_triangles; ++elem)
  {
    /* Compute local stiffness and load */
    get_local_stiffness(local_stiffness, &m, elem);
    get_local_load(local_load, &m, elem, fn_f);

#ifdef PRINT_DEBUG
    print_local_stiffness(local_stiffness);
    print_local_load(local_load);
#endif

    /* insert into global matrix and rhs */
    assemble_local2global_stiffness(local_stiffness, &mat, &m, elem);
    assemble_local2global_load(local_load, rhs, &m, elem);
  }

#ifdef PRINT_DEBUG
  printf("Matrix after assembly:\n");
  print_matrix(&mat);
#endif

  /* 4. Apply boundary conditions */
  apply_dbc(&mat, rhs, &m, fn_g);

#ifdef PRINT_DEBUG
  printf("Matrix after application of BCs:\n");
  print_matrix(&mat);
#endif

  /* 5. Solve the linear system */
  solve(&mat, u, rhs);

  /* 6. Evaluate error */
  errors[0] = l2_norm(u, &m, fn_u);
  errors[1] = inf_norm(u, &m, fn_u);

  /* free allocated resources */
  free(u);
  free(rhs);
  rhs = NULL;
  free_matrix(&mat);
  free_mesh(&m);
}

void print_local_stiffness(double local_stiffness[3][3])
{
  size_t i, j;
  for (i = 0; i < 3; ++i)
  {
    for (j = 0; j < 3; ++j)
    {
      printf("%5.2f ", local_stiffness[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_local_load(double local_load[3])
{
  size_t i;
  for (i = 0; i < 3; ++i)
  {
    printf("%5.2f ", local_load[i]);
  }
  printf("\n\n");
}

