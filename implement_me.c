/* -*- indent-tabs-mode: nil c-basic-offset: 2 -*- */
#include <stdio.h>
#include <string.h>

#include "exercise5.h"

#define err_exit(msg)                                                   \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

void* safe_malloc(size_t size)
{
  void *mem = malloc(size);
  if(NULL == mem) {
  fprintf(stderr, "malloc failed - exiting.\n");
  exit(EXIT_FAILURE);
  }
  memset(mem, 0, size);
  return mem;
}

/* add the vertex if it is not already there */
void add_neighbor(vertex_list **l, size_t v) {
  if(NULL == *l) {
    (*l) = safe_malloc(sizeof(vertex_list));
    (*l)->vertex = v;
    (*l)->next = NULL;
  }
  if((*l)->vertex != v) add_neighbor(&((*l)->next), v);
}

void print_list(vertex_list *l) {
  if(NULL == l) return;
  fprintf(stdout, " %lu", l->vertex);
  print_list(l->next);
}

void list_length(vertex_list *l) {
  if(NULL == l) return 0;
  return 1 + list_length(l->next);
}

enum boundary_t {
  SOUTH = 1,
  EAST  = 2,
  NORTH = 3,
  WEST  = 4
};

void get_mesh(mesh * m, size_t n)
{
  m->n_vertices = (n + 1) * (n + 1);
  m->n_triangles = 2 * n * n;
  m->coords      = safe_malloc(2 * m->n_vertices  * sizeof(double));
  m->t2v         = safe_malloc(3 * m->n_triangles * sizeof(size_t));
  m->id_v        = safe_malloc(m->n_vertices * sizeof(unsigned char));
  m->v_neighbors = safe_malloc(m->n_vertices * sizeof(vertex_list *));

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

      if(0 == iy) m->id_v[vertex_id] = (unsigned char) SOUTH;
      if(n == iy) m->id_v[vertex_id] = (unsigned char) NORTH;
      if(0 == ix) m->id_v[vertex_id] = (unsigned char) WEST;
      if(n == ix) m->id_v[vertex_id] = (unsigned char) EAST;
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
    add_neighbor(&m->v_neighbors[m->t2v[3 * tid + 0]], m->t2v[3 * tid + 1]);
    add_neighbor(&m->v_neighbors[m->t2v[3 * tid + 1]], m->t2v[3 * tid + 0]);
    add_neighbor(&m->v_neighbors[m->t2v[3 * tid + 1]], m->t2v[3 * tid + 2]);
    add_neighbor(&m->v_neighbors[m->t2v[3 * tid + 2]], m->t2v[3 * tid + 1]);
    add_neighbor(&m->v_neighbors[m->t2v[3 * tid + 2]], m->t2v[3 * tid + 0]);
    add_neighbor(&m->v_neighbors[m->t2v[3 * tid + 0]], m->t2v[3 * tid + 2]);
  }

  if(0) {
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
  if(NULL == m) mat = mat + 1 - 1;
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
