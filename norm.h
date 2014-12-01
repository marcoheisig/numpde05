#ifndef NORM_H
#define NORM_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mesh.h"

/*
 * Evaluates the discrete inf-norm at points given by the mesh
 */
double inf_norm(double const * u, mesh const * m,
                double (*u_exact)(double, double));

/*
 * Evaluates the discrete inf-norm at points given by the mesh
 * normalizes it by the number of points.
 */
double l2_norm(double const * u, mesh const * m,
               double (*u_exact)(double, double));

#ifdef __cplusplus
}
#endif

#endif
