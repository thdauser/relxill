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

#ifndef RELLP_EXTENDED_H_
#define RELLP_EXTENDED_H_

#include "Relbase.h"
#include "common.h"

/** dimensions of the emissivity table */
#define EXT_TABLE_NUM_SPINS 20 // number of spin values
#define EXT_TABLE_NUM_RADII 65 // number of height values
#define EXT_TABLE_NUM_ANGLES 18 // number of distance values
#define EXT_TABLE_NUM_RAD_BINS 99  // number of radial bin values
#define EXT_TABLE_NUM_ENERGY_BINS 20 // number of energy bin values
#define EXT_TABLE_FILENAME "rel_ext_table_disk_sph_v0.3.fits"

/** dimensions of the primary table */
#define PRIM_TABLE_NUM_SPINS 12 // number of spin values
#define PRIM_TABLE_NUM_RADII 65 // number of height values
#define PRIM_TABLE_NUM_ANGLES 18 // number of distance values
#define PRIM_TABLE_NUM_INCLINATIONS 20  // number of inclination bins
#define PRIM_TABLE_NUM_ENERGY_BINS 20 // number of energy bin values
#define PRIM_TABLE_FILENAME "rel_ext_table_obs_sph_v0.3.fits"

#define SLAB_NUMBER_RINGS 100

typedef struct {

  double* r_sph;
  double** theta;
  double* bin; // either radius or inclination
  double*** val; // either flux or lensing_factor
  double*** gshift;
  double** refl_frac;
  double** f_ad;
  double** f_inf;
  double** f_bh; // only primary table contains this value
  double**** energy_data_storage;
  // std::vector<std::vector<std::vector<std::vector<double>>>> energy_data_storage;

} extendedSourceData;

/** the extended source table structure */
typedef struct {
  float* a;
  int n_a;
  int n_r;
  int n_th;
  int n_bin;
  int n_storage;
  extendedSourceData** dat;
} extendedSourceTable;

typedef struct {
  int ind_a; // index of the spin
  double ifac_a; // interpolation factor between two spins
  int ind_r[2]; // two values - 0 is for lower spin of interpolation, 1 for higher spin
  double ifac_r[2];
  int ind_th[2];
  double ifac_th[2];
} InterpFactors;

typedef struct {
  int num_rings;
  VecD slab_grid_radii;
} slabParams;


// functions used globally:
void calc_emis_ring_source(emisProfile* emisProf, const relParam* param, int* status);

void calc_emis_slab_source(emisProfile* emisProf, const relParam* param, int* status);

double get_energyboost_ext_source_obs(const relParam* param, double gamma);

double get_energyboost_ext_source_disk(const relParam* param, double radius, double gamma);


// constructors and destructors:
extendedSourceTable* new_extended_source_table(int n_r, int n_th,
                                               int n_a, int n_bin, int n_storage,
                                               int* status);

extendedSourceData* new_extended_source_data(int* status, bool is_lensing_table);

void free_extended_source_table(extendedSourceTable* tab, bool is_lensing_table);

void free_extended_source_data(const extendedSourceData* dat, bool is_lensing_table);


#endif /* RELLP_EXTENDED_H_ */
