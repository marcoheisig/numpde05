/* -*- indent-tabs-mode: nil c-basic-offset: 2 -*- */
#pragma once
#include <stdlib.h>

/* I feel stupid for implementing this - 2014 and still builtin lists in C */
void* safe_malloc(size_t size);

struct node {
    size_t value;
    struct node *next;
};

typedef struct node node;

node *make_list();

void free_list(node *head);

/* add the vertex if not present, assert the list is ordered */
void add_ordered(node *head, size_t value);

void print_list1(node *head);

void print_list(node *head);

size_t list_length1(node *n, size_t length);

size_t list_length(node *head);
