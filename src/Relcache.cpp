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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "Relcache.h"
#include "Relbase.h"

/** probably best move to "utils" **/
inpar *get_inputvals_struct(double *ener, int n_ener, relParam *rel_par, int *status) {
  auto *inp = (inpar *) malloc(sizeof(inpar));
  CHECK_MALLOC_RET_STATUS(inp, status, nullptr)

  inp->ener = ener;
  inp->n_ener = n_ener;
  inp->rel_par = rel_par;

  return inp;
}

/** probably best move to "utils" **/
inpar *set_input_syspar(relParam *rel_par, int *status) {
  auto *inp = (inpar *) malloc(sizeof(inpar));
  CHECK_MALLOC_RET_STATUS(inp, status, nullptr)

  inp->ener = nullptr;
  inp->n_ener = 0;
  inp->rel_par = rel_par;

  return inp;
}

int comp_single_param_val(double val1, double val2) {
  if (fabs(val1 - val2) <= CACHE_LIMIT) {
    return 0;
  } else {
    return 1;
  }
}

static int comp_sys_param(const relParam *cpar, const relParam *par) {

  if (comp_single_param_val(par->a, cpar->a)) {
    return 1;
  }
  if (comp_single_param_val(par->emis1, cpar->emis1)) {
    return 1;
  }
  if (comp_single_param_val(par->emis2, cpar->emis2)) {
    return 1;
  }
  if (comp_single_param_val(par->gamma, cpar->gamma)) {
    return 1;
  }
  if (comp_single_param_val(par->height, cpar->height)) {
    return 1;
  }
  if (comp_single_param_val(par->htop, cpar->htop)) {
    return 1;
  }
  if (comp_single_param_val(par->incl, cpar->incl)) {
    return 1;
  }
  if (comp_single_param_val(par->beta, cpar->beta)) {
    return 1;
  }
  if (comp_single_param_val(par->rin, cpar->rin)) {
    return 1;
  }
  if (comp_single_param_val(par->rbr, cpar->rbr)) {
    return 1;
  }
  if (comp_single_param_val(par->rout, cpar->rout)) {
    return 1;
  }

  if (par->limb != cpar->limb) {
    return 1;
  }
  if (par->return_rad != cpar->return_rad) {
    return 1;
  }

  // for now, if we have correction factors, there is no caching of rel_param results possible
  // TODO: move correction factors and therefore the return rad emis profile outside of the system parameters
  if (par->rrad_corr_factors != nullptr || cpar->rrad_corr_factors!= nullptr ){
    return 1;
  }

  return 0;
}

int did_rel_param_change(const relParam *cpar, const relParam *par) {

  if (cpar == nullptr) {
    return 1;
  }

  // first check all system parameters
  if (comp_sys_param(cpar, par)) {
    return 1;
  }

  if (comp_single_param_val((double) par->emis_type, (double) cpar->emis_type)) {
    return 1;
  }
  if (comp_single_param_val((double) par->model_type, (double) cpar->model_type)) {
    return 1;
  }

  if (comp_single_param_val(par->z, cpar->z)) {
    return 1;
  }
  if (comp_single_param_val(par->lineE, cpar->lineE)) {
    return 1;
  }

  if (par->do_renorm_relline != cpar->do_renorm_relline) {
    return 1;
  }

  if (par->return_rad != cpar->return_rad) {
    return 1;
  }

  if (comp_single_param_val((double) par->ion_grad_type, (double) cpar->ion_grad_type)) return 1;

  /** also check if the number of zones changed **/
  if (par->num_zones != cpar->num_zones) {
    return 1;
  }

  return 0;
}

void set_cached_rel_param(relParam *par, relParam **ca_rel_param, int *status) {

  assert(ca_rel_param != nullptr);

  if ((*ca_rel_param) == nullptr) {
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
  (*ca_rel_param)->ion_grad_type = par->ion_grad_type;
  (*ca_rel_param)->num_zones = par->num_zones;


  (*ca_rel_param)->return_rad = par->return_rad;
  (*ca_rel_param)->rrad_corr_factors = par->rrad_corr_factors; // is not checked and therefore not used
}

void set_cached_xill_param(xillParam *par, xillParam **ca_xill_param, int *status) {

  if ((*ca_xill_param) == nullptr) {
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

  (*ca_xill_param)->iongrad_index = par->iongrad_index;

}

static int did_energy_grid_change(double *ener, int n_ener, relline_spec_multizone *ca) {
  int change = 0;

  if (ca == nullptr) {
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
  auto *new_node = new cnode;
  CHECK_MALLOC_RET_STATUS(new_node, status, nullptr)

  new_node->data = data;
  new_node->next = next;

  return new_node;
}

cnode *cli_prepend(cnode *head, cdata *data, int *status) {
  CHECK_STATUS_RET(*status, nullptr);
  cnode *new_node = cli_create(data, head, status);
  head = new_node;
  return head;
}


int cli_count_elements(cnode *head) {
  cnode *cursor = head;
  int nelem = 0;
  while (cursor != nullptr) {
    nelem++;
    cursor = cursor->next;
  }
  return nelem;
}

/* delete the linked list */
void cli_delete_list(cnode **pt_head) {
  cnode *cursor = *pt_head;

  cnode *next = nullptr;

  while (cursor != nullptr) {
    next = cursor->next;
    free_cnode(&cursor);
    cursor = next;
  }

  *pt_head = nullptr;
}

static cache_info *init_cache_info(cnode *node, int *status) {
  auto *ca = (cache_info *) malloc(sizeof(cache_info));
  CHECK_MALLOC_RET_STATUS(ca, status, nullptr)

  // set where the storage is pointing to, to the current node pointer
  ca->store = node;

  ca->read = nullptr;  // not used right now!!
  ca->store = nullptr;
  ca->relcache = 0;
  ca->syscache = 0;
  ca->xilcache = 0;

  return ca;
}

cnode *check_cache_syspar(cache_info *ca_info, inpar *inp, cnode *node) {

  if (comp_sys_param(node->data->par_rel, inp->rel_par) == 0) {
    // system parameters did not change in this iteration
    ca_info->syscache = 1;
    ca_info->read = nullptr;
    ca_info->store = node;
    return nullptr;

  } else {
    // parameters did change, let's try the next node
    return node->next;
  }
}

cnode *check_cache_relpar(cache_info *ca_info, inpar *inp, cnode *node) {

  if (did_rel_param_change(node->data->par_rel, inp->rel_par) == 0) {
    // system parameters did not change in this iteration
    // however, one last check if the energy grid did change
    if (did_energy_grid_change(inp->ener, inp->n_ener, node->data->relbase_spec)) {
      return node->next;

      // energy grid AND parameters did not change: found a MATCH
    } else {
      ca_info->relcache = 1;
      ca_info->read = nullptr;
      ca_info->store = node;
      return nullptr;
    }

  } else {
    // parameters did change, let's try the next node
    return node->next;
  }
}


static int get_cache_maxsize(){
  if (shouldOutfilesBeWritten()){
    return 1;
  } else {
    return CLI_NMAX;
  }

}

cache_info *cli_check_cache(cnode *head,
                            inpar *inp,
                            cnode *(*check_cache)(cache_info *, inpar *, cnode *),
                            int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  cache_info *ca_info = init_cache_info(head, status);
  CHECK_MALLOC_RET_STATUS(ca_info, status, nullptr)

  const int cache_maxsize = get_cache_maxsize();

  if (cache_maxsize==1){
    return ca_info;
  }

  int c = 0;
  cnode *cursor = head;
  cnode *next = nullptr;
  while (cursor != nullptr ) {

    // if cursor is not nullptr, we already have one element
    c++;

    // let's check the cache: return value can be nullptr for 2 conditions:
    // (1) found a match
    // (2) end of the list
    next = check_cache(ca_info, inp, cursor);

    // if we are above the maximal number of elements, delete the rest and break
    // (+) we need to set the cursor->next=nullptr
    if (next != nullptr && c >= cache_maxsize - 1) {
      if (is_debug_run()) {
        printf(" DEBUG: Cached reached its limiting size of %i\n", cache_maxsize);
      }
      cli_delete_list(&next);
      cursor->next = nullptr;
      assert(next == nullptr);
    }
    cursor = next;
  }

  return ca_info;
}

// prepend new node and set parameters
cnode *add_node_to_cache(cnode *head, relParam *relpar, xillParam *xillpar, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  cdata *data = init_cdata(status);

  assert(data != nullptr);
  if (relpar != nullptr) {
    set_cached_rel_param(relpar, &(data->par_rel), status);
  }

  if (xillpar != nullptr) {
    set_cached_xill_param(xillpar, &(data->par_xill), status);
  }

  // prepend (i.e., create new node, assign data, and set the new head)
  cnode *new_head = cli_prepend(head, data, status);

  assert(new_head != nullptr);

  CHECK_RELXILL_DEFAULT_ERROR(status);

  return new_head;
}

void add_relspec_to_cache(cnode **node, relParam *param, relline_spec_multizone *spec, int *status) {

  CHECK_STATUS_VOID(*status);

  cnode *old_head = *node;

  // prepend new node and set parameters
  cnode *new_head = add_node_to_cache(old_head, param, nullptr, status);

  // set the data
  new_head->data->relbase_spec = spec;

  *node = new_head;

  CHECK_RELXILL_DEFAULT_ERROR(status);
}

void set_cache_syspar(cnode **pt_head, relParam *param, RelSysPar *syspar, int *status) {

  CHECK_STATUS_VOID(*status);

  // prepend new node and set parameters
  cnode *new_head = add_node_to_cache(*pt_head, param, nullptr, status);

  // set the data
  new_head->data->relSysPar = syspar;

  *pt_head = new_head;

  CHECK_RELXILL_DEFAULT_ERROR(status);
}

/********* HELPER ROUTINES *********/


cdata *init_cdata(int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  auto *data = (cdata *) malloc(sizeof(cdata));
  CHECK_MALLOC_RET_STATUS(data, status, nullptr)

  data->par_rel = nullptr;
  data->par_xill = nullptr;
  data->relSysPar = nullptr;
  data->relbase_spec = nullptr;
  data->relxill_cache = nullptr;

  return data;
}

void free_relxill_cache(specCache *ca) {

  int ii;
  int m = 2;
  if (ca != nullptr) {
    if (ca->xill_spec != nullptr) {
      for (ii = 0; ii < ca->n_cache; ii++) {
        if (ca->xill_spec[ii] != nullptr) {
          free_xill_spec(ca->xill_spec[ii]);
        }
      }
      free(ca->xill_spec);
    }

    if (ca->fft_xill != nullptr) {
      free_fft_cache(ca->fft_xill, ca->n_cache, m);
    }

    if (ca->fft_rel != nullptr) {
      free_fft_cache(ca->fft_rel, ca->n_cache, m);
    }

    free_spectrum(ca->out_spec);

  }

  free(ca);

}

static void free_cdata(cdata **pt_data) {

  if (*pt_data != nullptr) {
    cdata *data = *pt_data;
    free(data->par_rel);
    free(data->par_xill);
    free_rel_spec(data->relbase_spec);
    free_relxill_cache(data->relxill_cache);
    free_relSysPar(data->relSysPar);
    free(data);
  }

}

void free_cnode(cnode **node) {

  if ((*node) != nullptr) {
    if ((*node)->data != nullptr) {
      free_cdata(&((*node)->data));
    }
    free(*node);
    *node = nullptr;
  }

}

int is_relbase_cached(cache_info *self) {
  return self->relcache;
}

int is_cached(cache_info *self) {
  return self->relcache + self->xilcache;
}
