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
  m->rank        = 0;
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
      else ++(m->rank);
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
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 0]], m->t2v[3 * tid + 1]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 1]], m->t2v[3 * tid + 0]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 1]], m->t2v[3 * tid + 2]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 2]], m->t2v[3 * tid + 1]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 2]], m->t2v[3 * tid + 0]);
    add_ordered(m->v_neighbors[m->t2v[3 * tid + 0]], m->t2v[3 * tid + 2]);
    /* Jep - each edge is added many times, but I don't care */
  }

  /* make an indexing scheme for interior nodes */
  m->id2row = safe_malloc(m->n_vertices * sizeof(size_t));
  m->row2id = safe_malloc(m->rank       * sizeof(size_t));
  size_t row = 0;
  for(size_t vid = 0; vid < m->n_vertices; ++vid) {
    if(m->id_v[vid] != INSIDE) continue;
    m->id2row[vid] = row;
    m->row2id[row] = vid;
    ++row;
  }

  if(1) {
    for(size_t vid = 0; vid < m->n_vertices; ++vid) {
      printf("neighbors of %lu:", vid);
      print_list(m->v_neighbors[vid]);
      printf("\n");
    }
  }

  /* output the mesh in gnuplot visualizable format */
  if(1) {
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
  }
}

void init_matrix(crs_matrix * mat, mesh const * m)
{
  mat->n_rows = m->rank;
  mat->n_cols = m->rank;

  size_t n_nonzero = 0;
  for(size_t row = 0; row < m->rank; ++row) {
    size_t n_connections = list_length(m->v_neighbors[m->row2id[row]]);
    n_nonzero += n_connections + 1;
  }
  mat->val    = safe_malloc(n_nonzero * sizeof(double));
  mat->colInd = safe_malloc(n_nonzero * sizeof(size_t));
  mat->rowPtr = safe_malloc(mat->n_rows * sizeof(size_t));

  //size_t index = 0;
  for(size_t row = 0; row < m->rank; ++row) {
    size_t n_connections = list_length(m->v_neighbors[m->row2id[row]]);
    for(size_t i = 0; i < n_connections; ++i) {
      
    }
  }
}

void get_local_stiffness(double local_stiffness[3][3], mesh const * m,
                         size_t element_id)
{
  err_exit("Not yet implemented!");
  /* Dummy operation to eliminate compiler errors */
  local_stiffness[0][0] = element_id + m->coords[0];
}

void get_local_load(double local_load[3], mesh const * m, size_t element_id,
                    double (*fn_f)(double, double))
{
  err_exit("Not yet implemented!");
  /* Dummy operation to eliminate compiler errors */
  local_load[0] = element_id + fn_f(m->coords[0], m->coords[1]);
}

void assemble_local2global_stiffness(double local_stiffness[3][3],
                                     crs_matrix * mat, mesh const * m,
                                     size_t element_id)
{
  err_exit("Not yet implemented!");
  /* Dummy operation to eliminate compiler errors */
  local_stiffness[0][0] = mat->val[0] + m->coords[0] + element_id;
}

void assemble_local2global_load(double local_load[3], double * rhs,
                                mesh const * m, size_t element_id)
{
  err_exit("Not yet implemented!");
  /* Dummy operation to eliminate compiler errors */
  local_load[0] = rhs[0] + m->coords[0] + element_id;
}

void apply_dbc(crs_matrix * mat, double * rhs, mesh const * m,
               double (fn_g)(unsigned char, double, double))
{
  err_exit("Not yet implemented!");
  /* Dummy operation to eliminate compiler errors */
  mat->val[0] = rhs[0] = fn_g(1, m->coords[0], m->coords[1]);
}

void solve(crs_matrix const * mat, double * u, double const * rhs)
{
  err_exit("Not yet implemented!");
  /* Dummy operation to eliminate compiler errors */
  u[0] = rhs[0] + mat->val[0];
}
