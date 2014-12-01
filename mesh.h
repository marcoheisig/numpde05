/* -*- indent-tabs-mode: nil c-basic-offset: 2 -*- */
#ifndef MESH_H
#define MESH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include "list.h"

  /*
 * A struct holding all relevant data structures to describe the triangulation.
 */
typedef struct {
  /* Interleaved node coordinates */
  double * coords;

  /* Vertex indices of each triangle in counter-clock-wise order */
  size_t * t2v;

  /* linked list of neighboring vertices of each vertex */
  node **v_neighbors;

  /* conversion between matrix row and mesh id */
  size_t *id2row;
  size_t *row2id;

  /* Vertex type id  (0: interior, 1: boundary y=0, 2: boundary x=1,
   * 3: boundary y=1, 4: boundary x=0) */
  unsigned char * id_v;

  /* Number of vertices */
  size_t n_vertices;

  /* rank = number of non-boundary vertices */
  size_t rank;

  /* Number of triangles */
  size_t n_triangles;
} mesh;

/*
 * Generates a triangulation of the unit square for a given n.
 */
void get_mesh(mesh * m, size_t n);

/*
 * Checks whether the given mesh is allocated.
 */
unsigned char is_allocated_mesh(mesh const * m);

/*
 * Prints the content of the mesh data structures to stdout.
 */
void print_mesh(mesh const * m);

/*
 * Frees all allocated resources of a mesh (but not the mesh object itself!).
 */
void free_mesh(mesh * m);

#ifdef __cplusplus
}
#endif

#endif
