#include <stdio.h>

#include "crs_matrix.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

unsigned char is_allocated_matrix(crs_matrix const * mat)
{
  return (mat != NULL && mat->n_rows > 0 && mat->n_cols > 0 &&
          mat->val != NULL && mat->rowPtr != NULL && mat->colInd != NULL);
}

void print_matrix(crs_matrix * mat)
{
  size_t i, j, ind;

  for (i = 0; i < mat->n_rows; ++i)
  {
    ind = mat->rowPtr[i];
    for (j = 0; j < mat->n_cols; ++j)
    {
      if (ind < mat->rowPtr[i+1] && mat->colInd[ind] == j)
        printf("%5.2f ", mat->val[ind++]);
      else
        printf("      ");
    }
    printf("\n");
  }
  printf("\n");
}

void free_matrix(crs_matrix * mat)
{
  free(mat->val);
  mat->val = NULL;
  free(mat->rowPtr);
  mat->rowPtr = NULL;
  free(mat->colInd);
  mat->colInd = NULL;
  mat->n_rows = 0;
  mat->n_cols = 0;
}
