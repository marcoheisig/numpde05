/* -*- indent-tabs-mode: nil c-basic-offset: 2 -*- */
#include <stdio.h>
#include <string.h>

#include "list.h"

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

node *make_list() {
  node *new_node = safe_malloc(sizeof(node));
  return new_node;
}

void free_list(node *l)
{
  if(NULL == l) return;
  free_list(l->next);
  free(l);
}

void add_ordered(node *head, size_t value)
{
  node *prev = head;
  node *current = head->next;

  while(current) {
    if(current->value == value) return;
    if(current->value > value) {
      node *new_node = safe_malloc(sizeof(node));
      new_node->value = value;
      new_node->next  = current;
      prev->next = new_node;
      return;
    }
    prev = current;
    current = current->next;
  }
  node *new_node = safe_malloc(sizeof(node));
  new_node->value = value;
  new_node->next  = current;
  prev->next = new_node;
}

void print_list1(node *l) {
  if(NULL == l) return;
  fprintf(stdout, " %lu", l->value);
  print_list1(l->next);
}

void print_list(node *head) {
  print_list1(head->next);
}

size_t list_length1(node *n, size_t length) {
  if(NULL == n) return length;
  return list_length1(n->next, length + 1);
}

size_t list_length(node *head) {
  return list_length1(head->next, 0);
}
