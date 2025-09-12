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
#define EXT_TABLE_NA 20 // number of spin values
#define EXT_TABLE_NHEIGHT 50 // number of height values
#define EXT_TABLE_NDISTANCE 17 // number of distance values
#define EXT_TABLE_NRAD 99  // number of radial bin values
#define EXT_TABLE_NENERGYBINS 20 // number of energy bin values
#define EXT_TABLE_FILENAME "rel_extended_source_sph.fits"

/** dimensions of the primary table */
#define PRIM_TABLE_NA 12 // number of spin values
#define PRIM_TABLE_NHEIGHT 50 // number of height values
#define PRIM_TABLE_NDISTANCE 17 // number of distance values
#define PRIM_TABLE_NINCL 20  // number of inclination bins
#define PRIM_TABLE_NENERGYBINS 20 // number of energy bin values
#define PRIM_TABLE_FILENAME "prim_extended_source_sph.fits"

typedef struct {

  double *height;
  // std::vector<std::vector<double>> x;
  double **x;
  double *bin; // was radius. Now either radius or inclination
  double ***val; // was flux. Now either flux or lensing_factor
  double ***gshift;
  double **refl_frac;
  double **f_ad;
  double **f_inf;
  double **f_bh; // now only primary table contains it
  double ****energy_data_storage;
  // std::vector<std::vector<std::vector<std::vector<double>>>> energy_data_storage;

} extendedSourceData;

/** the extended source table structure */
typedef struct {
  float *a;
  int n_a;
  int n_h;
  int n_x;
  int n_bin;
  int n_storage;
  extendedSourceData **dat;
} extendedSourceTable;

typedef struct {
  int ind_a;
  double ifac_a; // factor between two spins
  int ind_h[2]; // two values - 0 is for lower spin of interpolation, 1 is for higher spin
  double ifac_h[2];
  int ind_x[2];
  double ifac_x[2];
} InterpFactors;

// these I really need outside, all other below - arguable
void calc_emis_ring_source(emisProfile *emisProf, const relParam *param, int *status);
void calc_emis_disk_source(emisProfile *emisProf, const relParam *param, int *status);
double get_energyboost_ext_source_obs(const relParam *param, double gamma);
double get_energyboost_ext_source_disk(const relParam *param, double radius, double gamma);


// do I need to declare these two? I could just re-arrange code. Same applies to all below
void read_extendedSource_table(const char *filename, extendedSourceTable **inp_tab, int *status, bool is_lensing_table=false);
void calc_ring_emis(emisProfile *emisProf, double a, double h, double x, double gamma, extendedSourceTable *tab, int *status);
void calc_ring_fractions(emisProfile *emisProf, const relParam *param, extendedSourceTable *tab, int *status);


emisProfile *interpol_estable(extendedSourceTable *tab, InterpFactors ipol, double gamma, int *status);

// constructors and destructors
extendedSourceTable *new_extendedSourceTable(int n_h, int n_x,
                                             int n_a, int n_bin, int n_storage,
                                             int *status);

extendedSourceData *new_extendedSourceData(int *status, bool lensing);

void free_extendedSourceTable(extendedSourceTable *tab, bool lensing);

void free_extendedSourceData(extendedSourceData *dat, bool is_lensing_table);


#endif /* RELLP_EXTENDED_H_ */
