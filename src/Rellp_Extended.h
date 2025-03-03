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
#define ESTABLE_NA 20 // number of spin values
#define ESTABLE_NHEIGHT 50 // number of height values
#define ESTABLE_NDISTANCE 50 // number of distance values
#define ESTABLE_NRAD 99  // number of radial bin values
#define ESTABLE_NENERGYBINS 20 // number of energy bin values
#define ESTABLE_FILENAME "rel_extended_source_test.fits" // time to rename - this table is already final


typedef struct {

  double *height;
  double **r_ring;
  double *radius;
  double ***flux;
  double ***gshift;
  double **refl_frac;
  double **f_ad;
  double **f_inf;

  /* This will store the necessary table data (EnergyShift)
   * in order to sum over the 20 bins for a given Gamma exponent
   * which will then determine the gshift array
   * this array has dimensions NHEIGHT x NDISTANCE x NA x NRAD
   * x NENERGYBINS values corresponding to each radial bin
   */
  double ****energy_data_storage; // probably need to rename

} extendedSourceData;

/** the extended source table structure */
typedef struct {
  float *a;
  int n_a;
  int n_h;
  int n_rad;
  int n_storage;
  int n_r;
  extendedSourceData **dat;
} extendedSourceTable;


extendedSourceTable* get_es_table(int* status);

void printerror(int status);
void read_extendedSource_table(const char *filename, extendedSourceTable **inp_tab, int *status);
void get_emis_extendedSource(emisProfile *emisProf, const relParam *param, int *status);
void calc_extended_emis(emisProfile *emisProf, const relParam *param,
                        extendedSourceTable *tab, int *status);

// tbd
double calc_norm_factor_primary_spectrum_ext(double height, double a, double gamma, double f_inf);

emisProfile *interpol_estable(double a, double height, double r_ring,
                              double gamma, extendedSourceTable *tab, int *status);

void interpol_emisprofile_spin_height_distance(emisProfile *emis_profile_table,
                                               double ifac_a,
                                               const double *ifac_h,
                                               const double *ifac_r,
                                               const int *ind_h,
                                               const int *ind_r,
                                               extendedSourceData *const *dat,
                                               int nr);


// extendedSourceData functions
extendedSourceTable *new_extendedSourceTable(int n_h, int n_rad, int n_r,
                                             int n_a, int n_storage,
                                             int *status);

void free_extendedSourceData(extendedSourceData *dat);

extendedSourceData *new_extendedSourceData(int *status);

void free_extendedSourceTable(extendedSourceTable *tab);



#endif /* RELLP_EXTENDED_H_ */
