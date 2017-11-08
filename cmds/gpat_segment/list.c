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
#include "local_proto.h"
#include "list.h"

struct fifo* init_queue(unsigned index)
{
	struct fifo* fifo;
	struct node* node;

	fifo=malloc(sizeof(struct fifo));
	node=malloc(sizeof(struct node));
	node->index=index;
	node->next_node=NULL;
	fifo->first=node;
	fifo->last=node;
	fifo->length=1;

	return fifo;
}

int add_node(struct fifo* fifo, unsigned index)
{
	struct node* new_node;

	new_node=malloc(sizeof(struct node));
	new_node->next_node=NULL;
	new_node->index=index;
	new_node->count=1;
	fifo->last->next_node=new_node;
	fifo->last=new_node;
	fifo->length++;

	return 0;
}


int find_and_remove_node(struct fifo* fifo, unsigned index) {
	struct node* next;
	struct node* previous=NULL;

	next=fifo->first;

	while(next) {
		if(next->index!=index) {
			previous=next;
			next=next->next_node;
		} else {
			if(previous!=NULL)
				previous->next_node=next->next_node;
			else
				fifo->first=next->next_node;

			if(next==fifo->last) {
				fifo->last=previous;
				if(fifo->last)
					fifo->last->next_node=NULL;
			}
			free(next);
			fifo->length--;
			return fifo->length;
		}
	}
	return -1;
}

int get_value_at_position(struct fifo* fifo, unsigned position)
{
	int i=0;
	struct node* next;

	next=fifo->first;
	while(next) {
		if(i==position)
			return next->index;
		else {
			next=next->next_node;
			i++;
		}
	}

	return -1;
}

int in_queue(struct fifo* fifo, unsigned index)
{
	struct node* next;

	next=fifo->first;
	while(next) {
		if(next->index==index) {
			next->count++; /* this is side effect used to calculate the border*/
			return 1;
		}
		else
			next=next->next_node;
	}

	return 0;
}

int remove_node(struct fifo* fifo)
{
	struct node* old_node;

	old_node=fifo->first;
	fifo->first=old_node->next_node;
	free(old_node);
	fifo->length--;

	return 0;
}

unsigned pop_node(struct fifo* fifo)
{
	struct node* old_node;
	unsigned index;

	old_node=fifo->first;
	index=old_node->index;
	fifo->first=old_node->next_node;
	free(old_node);
	fifo->length--;

	return index;
}

int queue_length(struct fifo* fifo)
{
	if(fifo)
		return  fifo->length;
	else
		return 0;
}

int get_queue_length(struct fifo* fifo)
{
	int i=0;
	struct node* next=fifo->first;

	while(next) {
		next=next->next_node;
		i++;
	}

	return i;
}


int* convert_queue_to_table(struct fifo* fifo)
{
	int i=0;
	int* table;
	struct node* next;

	if(!fifo)
		return NULL;

	table=malloc(fifo->length*sizeof(int));
	next=fifo->first;

	get_queue_length(fifo);

	while(next) {
		table[i]=next->index;
		next=next->next_node;
		if(i>fifo->length)
			G_fatal_error("Cannot convert queue to table %d %d", i, fifo->length);
		i++;
	}

	return table;
}

int join_queues(struct fifo* a, struct fifo* b)
{
	a->last->next_node=b->first;
	a->last=b->last;
	a->length+=b->length;
	b->first=NULL;
	b->last=NULL;

	return 0;
}



int clear_queue(struct fifo* fifo)
{
	struct node* old_node;

	while(fifo->first) {
		old_node=fifo->first;
		fifo->first=old_node->next_node;
		free(old_node);
		fifo->length--;
	}

	return fifo->length;
}
