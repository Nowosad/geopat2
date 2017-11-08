#ifndef __LIST_H__
#define __LIST_H__

/****************************************************************************
 *
 * MODULE:	segmentation of motifels grid
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jacek Niesterowicz, Tomasz Stepinski
 * PURPOSE:	information retrieval using categorical maps:
 *		compares grid of histograms
 * COPYRIGHT:	(C) Space Informatics Lab, University of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

struct node{
	struct node* next_node;
	unsigned index;
	unsigned count;
};

struct fifo {
	struct node* first;
	struct node* last;
	int length;
};

struct fifo* init_queue(unsigned index);
int add_node(struct fifo* fifo, unsigned index);
int remove_node(struct fifo* fifo);
int in_queue(struct fifo* fifo, unsigned index);
int find_and_remove_node(struct fifo* fifo, unsigned index);
int get_value_at_position(struct fifo* fifo, unsigned position);
unsigned pop_node(struct fifo* fifo);
int queue_length(struct fifo* fifo);
int clear_queue(struct fifo* fifo);
int* convert_queue_to_table(struct fifo* fifo);
int join_queues(struct fifo* a, struct fifo* b);


/* add stack later */

#endif
