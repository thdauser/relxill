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
#ifndef RELCACHE_H_
#define RELCACHE_H_

#include "ModelDefinition.h"
#include "Xillspec.h"
#include <deque>
#include <utility>
#include "Xillspec.h"

extern "C" {
#include "common.h"
#include "relutility.h"
}

#define CACHE_LIMIT 1e-8

/****** TYPEDEF******/

#define CLI_NMAX 50
#define RELXILL_CACHE_SIZE 50

typedef struct cdata {

  xillParam *par_xill;
  relParam *par_rel;

  relline_spec_multizone *relbase_spec;
  RelSysPar *relSysPar;
  specCache *relxill_cache;
} cdata;

typedef struct cnode {
  cdata *data;
  struct cnode *next;
} cnode;

typedef struct cache_info {

  int relcache;
  int syscache;
  int xilcache;

  cnode *store;  // this is a pointer to the node, where the data is stored
  cnode *read;   // NOT NEEDED RIGHT NOW XXX


} cache_info;

typedef struct inpar {
  const relParam *rel_par;
  double *ener;
  int n_ener;
} inpar;

/** set the input parameters in one single structure **/
inpar *get_inputvals_struct(double *ener, int n_ener, const relParam *rel_par, int *status);
inpar *set_input_syspar(const relParam *rel_par, int *status);

bool are_values_different(double val1, double val2);

/** create a caching node **/
cnode *cli_create(cdata *data, cnode *next, int *status);

/** delete the full list (starting from head) **/
void cli_delete_list(cnode **head);
int cli_count_elements(cnode *head);

cnode *cli_prepend(cnode *head, cdata *data, int *status);
// void cli_traverse(cnode* head,callback f);

cdata *init_cdata(int *status);

// Routines to set the cached parameters
void add_relspec_to_cache(cnode **node, const relParam *param, relline_spec_multizone *spec, int *status);
void set_cache_syspar(cnode **node, const relParam *param, RelSysPar *syspar, int *status);
void set_cached_xill_param(xillParam *par, xillParam **ca_xill_param, int *status);

cnode *check_cache_syspar(cache_info *ca_info, inpar *input, cnode *node);
cnode *check_cache_relpar(cache_info *ca_info, inpar *input, cnode *node);

int did_rel_param_change(const relParam *cpar, const relParam *par);

cache_info *cli_check_cache(cnode *head,
                            inpar *inp,
                            cnode *(*check_cache)(cache_info *, inpar *, cnode *),
                            int *status);
//cache_info* cli_check_cache(cnode* head, input* input, cnode* *check_cache, int* status);

int is_relbase_cached(cache_info *self);
int is_xill_cached(cache_info *self);
int is_cached(cache_info *self);

void free_cnode(cnode **node);


class RelxillSpec : public Spectrum {

 public:
  RelxillSpec() : Spectrum(get_relxill_conv_energy_grid()->ener, get_relxill_conv_energy_grid()->nbins) {}

  explicit RelxillSpec(double *_flux) : RelxillSpec() {
    copy_flux(_flux);
  }

  // delete copy and move assignment constructor
  RelxillSpec &operator=(const RelxillSpec &other) = delete;

  RelxillSpec &operator=(const RelxillSpec &&other) = delete;

  // derived from the base class
  RelxillSpec(const RelxillSpec &inst) = default;

  RelxillSpec(RelxillSpec &&inst) = default;
};


class RelxillCacheElement {


 public:
  auto operator==(const RelxillCacheElement &_comp_cache) -> bool {
    auto parnames = m_model_definition.get_parnames();

    if (m_model_definition.get_model_name() != _comp_cache.m_model_definition.get_model_name()) {
      return false;
    }

     for(auto par=parnames.begin(); par!=parnames.cend(); ++par ) {
       if (are_values_different(m_model_definition.get_par(*par), _comp_cache.m_model_definition.get_par(*par))) {
         return false;
       }
     }
     return true;
  }

  bool model_params_identical(const ModelDefinition &_model_definition) {

    if (m_model_definition.get_model_name() != _model_definition.get_model_name()) {
      return false;
    }

    auto parnames = m_model_definition.get_parnames();
    for (auto par = parnames.begin(); par != parnames.cend(); ++par) {
      if (are_values_different(m_model_definition.get_par(*par), _model_definition.get_par(*par))) {
        return false;
      }
    }

    return true;
  }

  explicit RelxillCacheElement(ModelDefinition _model_definition, RelxillSpec _spec)
      : m_model_definition{_model_definition}, spec{(_spec)} {}


 private:
  ModelDefinition m_model_definition;
 public:
  const RelxillSpec spec;


};


class RelxillCache {

 public:


  /**
 * use this instance to always and securely access the database, without
 * needing to initialize it anywhere else
 * @return ModelDatabase
 */
  static RelxillCache &instance() {
    static auto *instance = new RelxillCache(RELXILL_CACHE_SIZE);
    return *instance;
  }


  void add(ModelDefinition _model_definition, RelxillSpec _spec) {

    if (m_cache.size() > max_size){
      m_cache.pop_back();
    }
    m_cache.emplace_front(_model_definition, _spec);

    if (is_debug_run() != 0) {
      printf(" Relxill Cache: added element, current size is %zu \n", m_cache.size());
    }
  }

  // empty the cache
  void empty() {
    while (!m_cache.empty()) {
      m_cache.pop_front();
    }
  }

  // get the size of the cache
  size_t size() {
    return m_cache.size();
  }

  [[nodiscard]] std::pair<bool, RelxillSpec> find_spec_pair(const ModelDefinition &_params) const {

    for (auto elem: m_cache) {
      if (elem.model_params_identical(_params)) {
        auto pair = std::make_pair(true, elem.spec);
        return pair;
      }
    }
    auto pair = std::make_pair(false, RelxillSpec());
    return pair;
  }

  [[nodiscard]] static RelxillSpec get_spec(const std::pair<bool, RelxillSpec> &cache_elem_pair) {
    return cache_elem_pair.second;
  }

  static bool is_cached(const std::pair<bool, RelxillSpec> &cache_elem_pair) {
    return cache_elem_pair.first;
  }


 private:
  explicit RelxillCache(size_t _max_size) : max_size{_max_size} {};

  std::deque<RelxillCacheElement> m_cache;
  size_t max_size;
};


#endif /* RELCACHE_H_ */
