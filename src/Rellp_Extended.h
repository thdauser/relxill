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

/** dimensions of the ES table */
#define EXT_TABLE_NA 20 // number of spin values
#define EXT_TABLE_NHEIGHT 50 // number of height values
#define EXT_TABLE_NDISTANCE 50 // number of distance values
#define EXT_TABLE_NRAD 99  // number of radial bin values
#define EXT_TABLE_NENERGYBINS 20 // number of energy bin values
#define EXT_TABLE_FILENAME "rel_extended_source_test.fits" // to rename, tests are gone

/** dimensions of the primary table */
#define PRIM_TABLE_NA 11 // number of spin values
#define PRIM_TABLE_NHEIGHT 50 // number of height values
#define PRIM_TABLE_NDISTANCE 50 // number of distance values
#define PRIM_TABLE_NINCL 20  // number of inclination bins
#define PRIM_TABLE_NENERGYBINS 20 // number of energy bin values
#define PRIM_TABLE_FILENAME "prim_extended_source.fits"

typedef struct {

  double *height;
  double **r_ring;
  double *bin; // was radius. Now either radius or inclination
  double ***val; // was flux. Now either flux or lensing
  double ***gshift;
  double **refl_frac;
  double **f_ad;
  double **f_inf;
  double **f_bh; // now only primary table contains it
  double ****energy_data_storage;

} extendedSourceData;

/** the extended source table structure */
typedef struct {
  float *a;
  int n_a;
  int n_h;
  int n_r;
  int n_bin;
  int n_storage;
  extendedSourceData **dat;
} extendedSourceTable;


extendedSourceTable* get_es_table(int* status);

void read_extendedSource_table(const char *filename, extendedSourceTable **inp_tab, int *status, bool lensing=false);
void get_emis_extendedSource(emisProfile *emisProf, const relParam *param, int *status);
void calc_extended_emis(emisProfile *emisProf, const relParam *param,
                        extendedSourceTable *tab, int *status);

emisProfile *interpol_estable(double a, double height, double r_ring,
                              double gamma, extendedSourceTable *tab, int *status);

// extendedSourceData functions
extendedSourceTable *new_extendedSourceTable(int n_h, int n_r,
                                             int n_a, int n_bin, int n_storage,
                                             int *status);

void free_extendedSourceData(extendedSourceData *dat, bool lensing);

extendedSourceData *new_extendedSourceData(int *status, bool lensing);

void free_extendedSourceTable(extendedSourceTable *tab, bool lensing);



#endif /* RELLP_EXTENDED_H_ */
