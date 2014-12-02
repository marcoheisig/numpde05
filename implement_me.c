/* -*- indent-tabs-mode: nil c-basic-offset: 2 -*- */
#include <stdio.h>
#include <string.h>

#include "exercise5.h"

#define err_exit(msg)                                                   \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

enum node_t {
  INSIDE = 0,
  SOUTH  = 1,
  EAST   = 2,
  NORTH  = 3,
  WEST   = 4
};

void get_mesh(mesh * m, size_t n)
{
  m->n_vertices = (n + 1) * (n + 1);
  m->n_triangles = 2 * n * n;
  m->coords      = safe_malloc(2 * m->n_vertices  * sizeof(double));
  m->t2v         = safe_malloc(3 * m->n_triangles * sizeof(size_t));
  m->id_v        = safe_malloc(m->n_vertices * sizeof(unsigned char));
  m->v_neighbors = safe_malloc(m->n_vertices * sizeof(node *));
  for(size_t vid = 0; vid < m->n_vertices; ++vid) {
    m->v_neighbors[vid] = make_list();
  }

  /* initialize vertices */
  for(size_t iy = 0; iy <= n; ++iy) {
    for(size_t ix = 0; ix <= n; ++ix) {
      double x, y;
      if(ix % 2 == 0) x = (double)ix / (double)n;
      else            x = (double)ix / (double)n + 1.0 / ((double)n * (double)n);
      if(iy % 2 == 0) y = (double)iy / (double)n;
      else            y = (double)iy / (double)n + 1.0 / ((double)n * (double)n);

      size_t vertex_id = iy * (n + 1) + ix;
      m->coords[2 * vertex_id + 0] = x;
      m->coords[2 * vertex_id + 1] = y;

      if     (0 == iy) m->id_v[vertex_id] = (unsigned char) SOUTH;
      else if(n == iy) m->id_v[vertex_id] = (unsigned char) NORTH;
      else if(0 == ix) m->id_v[vertex_id] = (unsigned char) WEST;
      else if(n == ix) m->id_v[vertex_id] = (unsigned char) EAST;
    }
  }

  /* initialize triangles by looping over all quadriliterals */
  for(size_t iy = 0; iy < n; ++iy) {
    for(size_t ix = 0; ix < n; ++ix) {
      // introduce local vertex ids within a quadriliteral
      size_t v00 = iy * (n + 1) + ix;
      size_t v10 = v00 + 1;
      size_t v01 = v00 + n + 1;
      size_t v11 = v01 + 1;

      size_t tid1 = 6 * (iy * n + ix);
      size_t tid2 = tid1 + 3;
      m->t2v[tid1 + 0] = v00;
      m->t2v[tid1 + 1] = v11;
      m->t2v[tid1 + 2] = v01;
      m->t2v[tid2 + 0] = v00;
      m->t2v[tid2 + 1] = v10;
      m->t2v[tid2 + 2] = v11;
    }
  }

  /* store connectivity information */
  for(size_t tid = 0; tid < m->n_triangles; ++tid) {
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 0]], m->t2v[3 * tid + 0]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 1]], m->t2v[3 * tid + 1]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 2]], m->t2v[3 * tid + 2]);

    add_ordered(m->v_neighbors[m->t2v[3 * tid + 0]], m->t2v[3 * tid + 1]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 1]], m->t2v[3 * tid + 0]);

    add_ordered(m->v_neighbors[m->t2v[3 * tid + 1]], m->t2v[3 * tid + 2]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 2]], m->t2v[3 * tid + 1]);

    add_ordered(m->v_neighbors[m->t2v[3 * tid + 2]], m->t2v[3 * tid + 0]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 0]], m->t2v[3 * tid + 2]);
    /* Jep - each edge is added many times, but I don't care */
  }

#ifdef PRINT_DEBUG
  for(size_t vid = 0; vid < m->n_vertices; ++vid) {
    printf("neighbors of %lu:", vid);
    print_list(m->v_neighbors[vid]);
    printf("\n");
  }
#endif

#ifdef PRINT_DEBUG
  /* output the mesh in gnuplot visualizable format */
  FILE *gnuplotfile = fopen("mesh.gnuplot", "w");
  if(!gnuplotfile) {
    fprintf(stderr, "could not open file for gnuplot visualization\n");
    return;
  }
  for(size_t tid = 0; tid < m->n_triangles; ++tid) {
    size_t v1 = m->t2v[3 * tid + 0];
    size_t v2 = m->t2v[3 * tid + 1];
    size_t v3 = m->t2v[3 * tid + 2];
    fprintf(gnuplotfile, "%g %g 0\n%g %g 0\n\n%g %g 0\n%g %g 0\n\n%g %g 0\n%g %g 0\n\n\n",
            m->coords[2 * v1 + 0], m->coords[2 * v1 + 1],
            m->coords[2 * v2 + 0], m->coords[2 * v2 + 1],
            m->coords[2 * v2 + 0], m->coords[2 * v2 + 1],
            m->coords[2 * v3 + 0], m->coords[2 * v3 + 1],
            m->coords[2 * v3 + 0], m->coords[2 * v3 + 1],
            m->coords[2 * v1 + 0], m->coords[2 * v1 + 1]);
  }
  fclose(gnuplotfile);
#endif
}

void init_matrix(crs_matrix * mat, mesh const * m)
{
  mat->n_rows = m->n_vertices;
  mat->n_cols = m->n_vertices;

  size_t n_values = 0;
  for(size_t row = 0; row < m->n_vertices; ++row) {
    node *neighbors = m->v_neighbors[row];
    node *current = neighbors->next; // skip head element - leaky abstraction
    while(current) {
      ++n_values;
      current = current->next;
    }
  }

  mat->val    = safe_malloc(n_values * sizeof(double));
  mat->colInd = safe_malloc(n_values * sizeof(size_t));
  mat->rowPtr = safe_malloc((1 + mat->n_rows) * sizeof(size_t));

  size_t index = 0;
  for(size_t row = 0; row < m->n_vertices; ++row) {
    mat->rowPtr[row] = index;
    node *neighbors = m->v_neighbors[row];
    node *current = neighbors->next;
    while(current) {
      mat->colInd[index] = current->value;
      ++index;
      current = current->next;
    }
  }
  mat->rowPtr[m->n_vertices] = index;
#ifdef PRINT_DEBUG
  printf("sparse matrix structure:\n");
  print_matrix(mat);
#endif
}

void get_local_stiffness(double local_stiffness[3][3], mesh const * m,
                         size_t element_id)
{
  double x1, x2, x3, y1, y2, y3;
  size_t v_id1 = m->t2v[3 * element_id + 0];
  size_t v_id2 = m->t2v[3 * element_id + 1];
  size_t v_id3 = m->t2v[3 * element_id + 2];
  x1 = m->coords[2 * v_id1 + 0];
  y1 = m->coords[2 * v_id1 + 1];
  x2 = m->coords[2 * v_id2 + 0];
  y2 = m->coords[2 * v_id2 + 1];
  x3 = m->coords[2 * v_id3 + 0];
  y3 = m->coords[2 * v_id3 + 1];

  double B11 = x2 - x1;
  double B12 = x3 - x1;
  double B21 = y2 - y1;
  double B22 = y3 - y1;
  double detB = fabs(B11 * B22 - B21 * B12);

#ifdef PRINT_DEBUG
  printf("generating local stiffness matrix for triangle %zu\n", element_id);
  printf("vertices: %zu %zu %zu\n", v_id1, v_id2, v_id3);
  printf("B :=\n| %g-%g %g-%g |\n| %g-%g %g-%g |\n",
         x2, x1, x3, x1, y2, y1, y3, y1);
  printf("B :=\n| %g %g |\n| %g %g |\n",
         B11, B12, B21, B22);
  printf("detB: %g\n\n", detB);
#endif

  double gamma1 = 1.0 / detB * ((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1));
  double gamma2 = 1.0 / detB * ((x2-x1)*(x3-x1) + (y2-y1)*(y3-y1));
  double gamma3 = 1.0 / detB * ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

  // contribution from S1
  local_stiffness[0][0] +=  0.5 * gamma1;
  local_stiffness[0][1] += -0.5 * gamma1;
  local_stiffness[1][0] += -0.5 * gamma1;
  local_stiffness[1][1] +=  0.5 * gamma1;

  // contribution from S2
  local_stiffness[0][0] +=  1.0 * gamma2;
  local_stiffness[0][1] += -0.5 * gamma2;
  local_stiffness[0][2] += -0.5 * gamma2;
  local_stiffness[1][0] += -0.5 * gamma2;
  local_stiffness[1][0] += -0.5 * gamma2;
  local_stiffness[2][1] +=  0.5 * gamma2;
  local_stiffness[1][2] +=  0.5 * gamma2;

  // contribution from S3
  local_stiffness[0][0] +=  0.5 * gamma3;
  local_stiffness[0][2] += -0.5 * gamma3;
  local_stiffness[2][0] += -0.5 * gamma3;
  local_stiffness[2][2] +=  0.5 * gamma3;
}

void get_local_load(double local_load[3], mesh const * m, size_t element_id,
                    double (*fn_f)(double, double))
{
  if(m == NULL) return;
  if(element_id == 0) if(fn_f == NULL) return;

  local_load[0] = 0.0;
  local_load[1] = 0.0;
  local_load[2] = 0.0;
}

void assemble_local2global_stiffness(double local_stiffness[3][3],
                                     crs_matrix * mat, mesh const * m,
                                     size_t element_id)
{
  for(size_t i = 0; i < 3; ++i) {
    size_t row = m->t2v[3*element_id + i];
    for(size_t j = 0; j < 3; ++j) {
      size_t col = m->t2v[3*element_id + j];
      for(size_t crs_i = mat->rowPtr[row]; crs_i < mat->rowPtr[row+1]; ++crs_i) {
        if(mat->colInd[crs_i] != col) continue;
        mat->val[crs_i] += local_stiffness[i][j];
        break;
      }
    }
  }
}

void assemble_local2global_load(double local_load[3], double * rhs,
                                mesh const * m, size_t element_id)
{
  for(size_t i = 0; i < 3; ++i) {
    size_t row = m->t2v[3*element_id + i];
    rhs[row] += local_load[i];
  }
}

void apply_dbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_g)(unsigned char, double, double))
{
  for(size_t i = 0; i < mat->n_rows; ++i) {
    if(m->id_v[i] != INSIDE) {
      double x = m->coords[i + 0];
      double y = m->coords[i + 1];
      for(size_t crsi = mat->rowPtr[i]; crsi < mat->rowPtr[i+1]; ++crsi) {
        size_t j = mat->colInd[crsi];
        if(i == j) {
          mat->val[crsi] = 1.0;
          rhs[i] = fn_g(m->id_v[i], x, y);
        } else {
          rhs[i] -= mat->val[crsi] * fn_g(m->id_v[i], x, y);
          mat->val[crsi] = 0.0;
        }
      }
    }
  }
}

void solve(crs_matrix const * mat, double * u, double const * rhs)
{
  err_exit("Not yet implemented!");
  /* Dummy operation to eliminate compiler errors */
  u[0] = rhs[0] + mat->val[0];
}
