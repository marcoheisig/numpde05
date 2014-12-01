#ifndef FEM_H
#define FEM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include "mesh.h"
#include "crs_matrix.h"

/*
 * Solver routine which carries out the FEM algorithm on a triangular grid:
 * - Preprocessing (Triangulation)
 * - Assembly of the LSE
 * - Application of BCs
 * - Solving of the LSE
 *
 * Parameters:
 * - n is the number of rectangle cells per spatial dimension, into which the
 *   unit square is split.
 * - errors is a 2-element double precision array, which will hold afterwards
 *   the deviation from the exact solution measured in the L2-norm and inf-norm.
 * - fn_f is a function pointer to the rhs function.
 * - fn_g is a function pointer to the Dirichlet boundary value function.
 * - fn_u is a function pointer to the analytical solution.
 */
void fem(size_t n, double errors[2], double (*fn_f)(double, double),
         double (*fn_g)(unsigned char, double, double),
         double (*fn_u)(double, double));

/*
 * Computes the local stiffness matrix for the given element number
 */
void get_local_stiffness(double local_stiffness[3][3], mesh const * m,
                         size_t element_id);

/*
 * Computes the local load vector for the given element number, using the
 * given function for the right hand side.
 */
void get_local_load(double local_load[3], mesh const * m, size_t element_id,
                    double (*fn_f)(double, double));

/*
 * Assembles the contributions of a local into the global stiffness matrix.
 */
void assemble_local2global_stiffness(double local_stiffness[3][3],
                                     crs_matrix * mat, mesh const * m,
                                     size_t element_id);

/*
 * Assembles the contributions of a local to the global load vector.
 */
void assemble_local2global_load(double local_load[3], double * rhs,
                                mesh const * m, size_t element_id);

/*
 * Applies the Dirichlet boundary conditions to the linear system.
 */
void apply_dbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_g)(unsigned char i, double, double));

/*
 * Solves a linear system Au=f.
 */
void solve(crs_matrix const * mat, double * u, double const * rhs);

/*
 * Prints a local stiffness matrix
 */
void print_local_stiffness(double local_stiffness[3][3]);

/*
 * Prints a local load vector
 */
void print_local_load(double local_load[3]);

#ifdef __cplusplus
}
#endif

#endif
