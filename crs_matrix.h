#ifndef CRS_MATRIX_H
#define CRS_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include "mesh.h"

/*
 * A struct holding all relevant data structures for a CRS-matrix.
 */
typedef struct {
  /* Values of matrix entries */
  double * val;

  /* Offset for each row in the storage arrays val and colInd */
  size_t * rowPtr;

  /* Column index for each stored value */
  size_t * colInd;

  /* Number of rows */
  size_t n_rows;

  /* Number of columns */
  size_t n_cols;
} crs_matrix;

/*
 * Allocates and initializes the CRS-matrix according to the connectivity given
 * in the mesh object. All entries are initialized to zero.
 */
void init_matrix(crs_matrix * mat, mesh const * m);

/*
 * Checks whether the given matrix is allocated.
 */
unsigned char is_allocated_matrix(crs_matrix const * mat);

/*
 * Prints the entire matrix to stdout
 */
void print_matrix(crs_matrix * mat);

/*
 * Frees all allocated resources of a matrix (but not the matrix object itself!).
 */
void free_matrix(crs_matrix * mat);


#ifdef __cplusplus
}
#endif

#endif
