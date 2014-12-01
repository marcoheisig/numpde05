/* -*- indent-tabs-mode: nil c-basic-offset: 2 -*- */
#include <stdio.h>

#include "mesh.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

unsigned char is_allocated_mesh(mesh const * m)
{
  return (m != NULL && m->n_vertices > 0 && m->n_triangles > 0 &&
          m->coords != NULL && m->t2v != NULL && m->id_v != NULL);
}

void print_mesh(mesh const * m)
{
  size_t i;

  if (!is_allocated_mesh(m)) err_exit("Mesh is not allocated");

  printf("Vertex coordinates: ");
  for (i = 0; i < m->n_vertices; ++i)
    printf("(%5.2f, %5.2f) ", m->coords[2*i], m->coords[2*i+1]);
  printf("\n\n");

  printf("Vertices per triangle: ");
  for (i = 0; i < m->n_triangles; ++i)
    printf("(%d: %d, %d, %d) ", (int)i, (int)m->t2v[3*i], (int)m->t2v[3*i+1], (int)m->t2v[3*i+2]);
  printf("\n\n");

  printf("Vertex ids: ");
  for (i = 0; i < m->n_vertices; ++i)
    printf("(%d: %d) ", (int)i, (int)m->id_v[i]);
  printf("\n\n");
}

void free_mesh(mesh * m)
{
  free(m->coords);
  m->coords = NULL;
  free(m->t2v);
  m->t2v = NULL;
  free(m->id_v);
  m->id_v = NULL;
  m->n_vertices = 0;
  m->n_triangles = 0;

  /* free the linked list of neighbors */
  for(size_t i = 0; i < m->n_vertices; ++i) {
    free_list(m->v_neighbors[i]);
  }
  free(m->v_neighbors);
}
