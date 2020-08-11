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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "relcache.h"

/** probably best move to "utils" **/
inpar *set_input(double *ener, int n_ener, relParam *rel_par, xillParam *xill_par, int *status) {
  inpar *inp = (inpar *) malloc(sizeof(inpar));
  CHECK_MALLOC_RET_STATUS(inp, status, NULL)

  inp->ener = ener;
  inp->n_ener = n_ener;
  inp->rel_par = rel_par;
  inp->xil_par = xill_par;

  return inp;
}

/** probably best move to "utils" **/
inpar *set_input_syspar(relParam *rel_par, int *status) {
  inpar *inp = (inpar *) malloc(sizeof(inpar));
  CHECK_MALLOC_RET_STATUS(inp, status, NULL)

  inp->ener = NULL;
  inp->n_ener = 0;
  inp->rel_par = rel_par;
  inp->xil_par = NULL;

  return inp;
}

int comp_single_param_val(double val1, double val2) {
  if (fabs(val1 - val2) <= CACHE_LIMIT) {
    return 0;
  } else {
    return 1;
  }
}

static int comp_sys_param(relParam *cpar, relParam *par) {

  if (comp_single_param_val(par->a, cpar->a)) return 1;
  if (comp_single_param_val(par->emis1, cpar->emis1)) return 1;
  if (comp_single_param_val(par->emis2, cpar->emis2)) return 1;
  if (comp_single_param_val(par->gamma, cpar->gamma)) return 1;
  if (comp_single_param_val(par->height, cpar->height)) return 1;
  if (comp_single_param_val(par->htop, cpar->htop)) return 1;
  if (comp_single_param_val(par->incl, cpar->incl)) return 1;
  if (comp_single_param_val(par->beta, cpar->beta)) return 1;
  if (comp_single_param_val(par->rin, cpar->rin)) return 1;
  if (comp_single_param_val(par->rbr, cpar->rbr)) return 1;
  if (comp_single_param_val(par->rout, cpar->rout)) return 1;
  if (par->limb != cpar->limb) return 1;

  return 0;
}

int comp_rel_param(relParam *cpar, relParam *par) {

  if (cpar == NULL) return 1;

  // first check all system parameters
  if (comp_sys_param(cpar, par)) return 1;

  if (comp_single_param_val((double) par->emis_type, (double) cpar->emis_type)) return 1;
  if (comp_single_param_val((double) par->model_type, (double) cpar->model_type)) return 1;

  if (comp_single_param_val(par->z, cpar->z)) return 1;
  if (comp_single_param_val(par->lineE, cpar->lineE)) return 1;

  if (par->do_renorm_relline != cpar->do_renorm_relline) return 1;

  /** also check if the number of zones changed **/
  if (par->num_zones != cpar->num_zones) return 1;

  return 0;
}

void set_cached_rel_param(relParam *par, relParam **ca_rel_param, int *status) {

  assert(ca_rel_param != NULL);

  if ((*ca_rel_param) == NULL) {
    (*ca_rel_param) = (relParam *) malloc(sizeof(relParam));
    CHECK_MALLOC_VOID_STATUS((*ca_rel_param), status)
  }

  (*ca_rel_param)->a = par->a;
  (*ca_rel_param)->emis1 = par->emis1;
  (*ca_rel_param)->emis2 = par->emis2;
  (*ca_rel_param)->gamma = par->gamma;
  (*ca_rel_param)->height = par->height;
  (*ca_rel_param)->htop = par->htop;
  (*ca_rel_param)->incl = par->incl;
  (*ca_rel_param)->beta = par->beta;

  (*ca_rel_param)->z = par->z;
  (*ca_rel_param)->limb = par->limb;
  (*ca_rel_param)->lineE = par->lineE;

  (*ca_rel_param)->emis_type = par->emis_type;
  (*ca_rel_param)->model_type = par->model_type;

  (*ca_rel_param)->rbr = par->rbr;
  (*ca_rel_param)->rin = par->rin;
  (*ca_rel_param)->rout = par->rout;

  (*ca_rel_param)->do_renorm_relline = par->do_renorm_relline;
  (*ca_rel_param)->num_zones = par->num_zones;

}

void set_cached_xill_param(xillParam *par, xillParam **ca_xill_param, int *status) {

  if ((*ca_xill_param) == NULL) {
    (*ca_xill_param) = (xillParam *) malloc(sizeof(xillParam));
    CHECK_MALLOC_VOID_STATUS((*ca_xill_param), status)
  }

  (*ca_xill_param)->afe = par->afe;
  (*ca_xill_param)->dens = par->dens;
  (*ca_xill_param)->ect = par->ect;
  (*ca_xill_param)->gam = par->gam;
  (*ca_xill_param)->lxi = par->lxi;
  (*ca_xill_param)->kTbb = par->kTbb;
  (*ca_xill_param)->frac_pl_bb = par->frac_pl_bb;

  (*ca_xill_param)->z = par->z;

  (*ca_xill_param)->prim_type = par->prim_type;
  (*ca_xill_param)->model_type = par->model_type;

  (*ca_xill_param)->ion_grad_index = par->ion_grad_index;
  (*ca_xill_param)->ion_grad_type = par->ion_grad_type;

}

static int did_energy_grid_change(double *ener, int n_ener, rel_spec *ca) {
  int change = 0;

  if (ca == NULL) {
    return change;
  }

  if (n_ener != ca->n_ener) {
    return 1;
  }

  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    if (fabs(ca->ener[ii] - ener[ii]) > 1e-4) {
      return 1;
    }
  }

  return change;
}

cnode *cli_create(cdata *data, cnode *next, int *status) {
  cnode *new_node = (cnode *) malloc(sizeof(cnode));
  CHECK_MALLOC_RET_STATUS(new_node, status, NULL)

  new_node->data = data;
  new_node->next = next;

  return new_node;
}

cnode *cli_prepend(cnode *head, cdata *data, int *status) {
  CHECK_STATUS_RET(*status, NULL);
  cnode *new_node = cli_create(data, head, status);
  head = new_node;
  return head;
}


/* traverse the linked list */
/*void cli_traverse(cnode* head,callback f){
    cnode* cursor = head;
    while(cursor != NULL) {
        f(cursor);
        cursor = cursor->next;
    }
}*/

int cli_count_elements(cnode *head) {
  cnode *cursor = head;
  int nelem = 0;
  while (cursor != NULL) {
    nelem++;
    cursor = cursor->next;
  }
  return nelem;
}

/* delete the linked list */
void cli_delete_list(cnode **pt_head) {
  cnode *cursor = *pt_head;

  cnode *next = NULL;

  while (cursor != NULL) {
    next = cursor->next;
    free_cnode(&cursor);
    cursor = next;
  }

  *pt_head = NULL;
}

static cache_info *init_cache_info(cnode *node, int *status) {
  cache_info *ca = (cache_info *) malloc(sizeof(cache_info));
  CHECK_MALLOC_RET_STATUS(ca, status, NULL)

  // set where the storage is pointing to, to the current node pointer
  ca->store = node;

  ca->read = NULL;  // not used right now!!
  ca->store = NULL;
  ca->relcache = 0;
  ca->syscache = 0;
  ca->xilcache = 0;

  return ca;
}

cnode *check_cache_syspar(cache_info *ca_info, inpar *inp, cnode *node) {

  if (comp_sys_param(node->data->par_rel, inp->rel_par) == 0) {
    // system parameters did not change in this iteration
    ca_info->syscache = 1;
    ca_info->read = NULL;
    ca_info->store = node;
    return NULL;

  } else {
    // parameters did change, let's try the next node
    return node->next;
  }
}

cnode *check_cache_relpar(cache_info *ca_info, inpar *inp, cnode *node) {

  if (comp_rel_param(node->data->par_rel, inp->rel_par) == 0) {
    // system parameters did not change in this iteration
    // however, one last check if the energy grid did change
    if (did_energy_grid_change(inp->ener, inp->n_ener, node->data->relbase_spec)) {
      return node->next;

      // energy grid AND parameters did not change: found a MATCH
    } else {
      ca_info->relcache = 1;
      ca_info->read = NULL;
      ca_info->store = node;
      return NULL;
    }

  } else {
    // parameters did change, let's try the next node
    return node->next;
  }
}

cache_info *cli_check_cache(cnode *head,
                            inpar *inp,
                            cnode *(*check_cache)(cache_info *, inpar *, cnode *),
                            int *status) {

  CHECK_STATUS_RET(*status, NULL);

  cache_info *ca_info = init_cache_info(head, status);
  CHECK_MALLOC_RET_STATUS(ca_info, status, NULL)

  int c = 0;
  cnode *cursor = head;
  cnode *next = NULL;
  while (cursor != NULL) {

    // if cursor is not NULL, we already have one element
    c++;

    // let's check the cache: return value can be NULL for 2 conditions:
    // (1) found a match
    // (2) end of the list
    next = check_cache(ca_info, inp, cursor);

    // if we are above the maximal number of elements, delete the rest and break
    // (+) we need to set the cursor->next=NULL
    if (next != NULL && c >= CLI_NMAX - 1) {
      if (is_debug_run()) {
        printf(" DEBUG: Cached reached its limiting size of %i\n", CLI_NMAX);
      }
      cli_delete_list(&next);
      cursor->next = NULL;
      assert(next == NULL);
    }
    cursor = next;
  }

  return ca_info;
}

// prepend new node and set parameters
cnode *add_node_to_cache(cnode *head, relParam *relpar, xillParam *xillpar, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  cdata *data = init_cdata(status);

  assert(data != NULL);
  if (relpar != NULL) {
    set_cached_rel_param(relpar, &(data->par_rel), status);
  }

  if (xillpar != NULL) {
    set_cached_xill_param(xillpar, &(data->par_xill), status);
  }

  // prepend (i.e., create new node, assign data, and set the new head)
  cnode *new_head = cli_prepend(head, data, status);

  assert(new_head != NULL);

  CHECK_RELXILL_DEFAULT_ERROR(status);

  return new_head;
}

void set_cache_relbase(cnode **pt_head, relParam *param, rel_spec *spec, int *status) {

  CHECK_STATUS_VOID(*status);

  cnode *old_head = *pt_head;

  // prepend new node and set parameters
  cnode *new_head = add_node_to_cache(old_head, param, NULL, status);

  // set the data
  new_head->data->relbase_spec = spec;

  *pt_head = new_head;

  CHECK_RELXILL_DEFAULT_ERROR(status);
}

void set_cache_syspar(cnode **pt_head, relParam *param, relSysPar *syspar, int *status) {

  CHECK_STATUS_VOID(*status);

  // prepend new node and set parameters
  cnode *new_head = add_node_to_cache(*pt_head, param, NULL, status);

  // set the data
  new_head->data->relSysPar = syspar;

  *pt_head = new_head;

  CHECK_RELXILL_DEFAULT_ERROR(status);
}

/********* HELPER ROUTINES *********/


cdata *init_cdata(int *status) {

  CHECK_STATUS_RET(*status, NULL);

  cdata *data = (cdata *) malloc(sizeof(cdata));
  CHECK_MALLOC_RET_STATUS(data, status, NULL)

  data->par_rel = NULL;
  data->par_xill = NULL;
  data->relSysPar = NULL;
  data->relbase_spec = NULL;
  data->relxill_cache = NULL;

  return data;
}

void free_relxill_cache(specCache *ca) {

  int ii;
  int m = 2;
  if (ca != NULL) {
    if (ca->xill_spec != NULL) {
      for (ii = 0; ii < ca->n_cache; ii++) {
        if (ca->xill_spec[ii] != NULL) {
          free_xill_spec(ca->xill_spec[ii]);
        }
      }
      free(ca->xill_spec);
    }

    if (ca->fft_xill != NULL) {
      free_fft_cache(ca->fft_xill, ca->n_cache, m);
    }

    if (ca->fft_rel != NULL) {
      free_fft_cache(ca->fft_rel, ca->n_cache, m);
    }

    free_out_spec(ca->out_spec);

  }

  free(ca);

}

static void free_cdata(cdata **pt_data) {

  if (*pt_data != NULL) {
    cdata *data = *pt_data;
    free(data->par_rel);
    free(data->par_xill);
    free_rel_spec(data->relbase_spec);
    free_relxill_cache(data->relxill_cache);
    free_relSysPar(data->relSysPar);
  }

}

void free_cnode(cnode **node) {

  if ((*node) != NULL) {
    if ((*node)->data != NULL) {
      free_cdata(&((*node)->data));
    }
    free(*node);
    *node = NULL;
  }

}

int is_relbase_cached(cache_info *self) {
  return self->relcache;
}

int is_xill_cached(cache_info *self) {
  return self->xilcache;
}

int is_cached(cache_info *self) {
  return self->relcache + self->xilcache;
}
