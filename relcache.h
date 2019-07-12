/*
   This file is part of the RELXILL model code.

   RELXILL is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   RELXILL is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2019 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELCACHE_H_
#define RELCACHE_H_

#include "common.h"
#include "relutility.h"

#define CACHE_LIMIT 1e-8
//const double cache_limit = 1e-8;


/****** TYPEDEF******/

#define CLI_NMAX 10


typedef struct cdata{

	xillParam* par_xill;
	relParam* par_rel;

	rel_spec* relbase_spec;
	relSysPar* relSysPar;
	specCache* relxill_cache;
} cdata;

typedef struct cnode{
    cdata* data;
    struct cnode* next;
} cnode;



typedef struct cache_info{

	int relcache;
	int syscache;
	int xilcache;

	cnode* store;  // this is a pointer to the node, where the data is stored
	cnode* read;   // NOT NEEDED RIGHT NOW XXX


} cache_info;




typedef struct inpar{

	relParam* rel_par;
	xillParam* xil_par;
	double* ener;
	int n_ener;
} inpar;

//typedef void (*callback)(cache_info* ca_info, input* input, cnode* node);

/** set the input parameters in one single structure **/
inpar* set_input(double* ener,int n_ener,relParam* rel_par,xillParam* xill_par, int* status);
inpar* set_input_syspar(relParam* rel_par, int* status);

int comp_single_param_val(double val1, double val2);

/** create a caching node **/
cnode* cli_create(cdata* data,cnode* next, int* status);

/** delete the full list (starting from head) **/
void cli_delete_list(cnode** head);
int cli_count_elements(cnode* head);

cnode* cli_prepend(cnode* head,cdata* data, int* status);
// void cli_traverse(cnode* head,callback f);

cdata* init_cdata(int* status);

// Routines to set the cached parameters
void set_cache_relbase(cnode** node, relParam* param, rel_spec* spec, int* status);
void set_cache_syspar(cnode** node, relParam* param, relSysPar* syspar, int* status);
void set_cached_xill_param(xillParam* par, xillParam** ca_xill_param, int* status);

cnode* check_cache_syspar(cache_info* ca_info, inpar* input, cnode* node);
cnode* check_cache_relpar(cache_info* ca_info, inpar* input, cnode* node);

int comp_rel_param(relParam* cpar, relParam* par);

cache_info* cli_check_cache(cnode* head, inpar* inp, cnode* (*check_cache) (cache_info*, inpar*, cnode*) , int* status);
//cache_info* cli_check_cache(cnode* head, input* input, cnode* *check_cache, int* status);

int is_relbase_cached(cache_info* self);
int is_xill_cached(cache_info* self);
int is_cached(cache_info* self);

void free_cnode(cnode** node);


#endif /* RELCACHE_H_ */
