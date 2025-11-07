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
#define EXT_TABLE_NUM_RADII 50 // number of height values
#define EXT_TABLE_NUM_ANGLES 17 // number of distance values
#define EXT_TABLE_NUM_RAD_BINS 99  // number of radial bin values
#define EXT_TABLE_NUM_ENERGY_BINS 20 // number of energy bin values
#define EXT_TABLE_FILENAME "rel_ext_table_disk_sph_v0.1.fits"

/** dimensions of the primary table */
#define PRIM_TABLE_NUM_SPINS 12 // number of spin values
#define PRIM_TABLE_NUM_RADII 50 // number of height values
#define PRIM_TABLE_NUM_ANGLES 17 // number of distance values
#define PRIM_TABLE_NUM_INCLINATIONS 20  // number of inclination bins
#define PRIM_TABLE_NUM_ENERGY_BINS 20 // number of energy bin values
#define PRIM_TABLE_FILENAME "rel_ext_table_obs_sph_v0.1.fits"

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

// functions used only in Rellp_Extended:
extendedSourceTable* get_es_table(int* status);
extendedSourceTable* get_prim_table(int* status);

static void get_spin_rows_extended_table(const int nrows, float** val, fitsfile* fptr, int* status);
static extendedSourceData* load_single_extended_source_data(
          fitsfile* fptr, const int n_h, const int n_x, const int n_bin, const int n_storage,
          int* status, const double spin, const bool is_lensing_table);
void read_extended_source_table(const char* filename, extendedSourceTable** inp_tab, int* status,
                                const bool is_lensing_table=false);

emisProfile* interpol_extended_table(extendedSourceTable* tab, InterpFactors ipol,
                                     const double gamma, double* radial_grid, int* status);

double get_interpolated_value_a_h_x(const double ifac_a, const double ifac_h[2], const double ifac_x[2],
                                    const double vals_lo00, const double vals_lo01, const double vals_lo10, const double vals_lo11,
                                    const double vals_hi00, const double vals_hi01, const double vals_hi10, const double vals_hi11);
void get_ipol_factors_a_h_x(double a, double h, double x,
                            const extendedSourceTable* tab, InterpFactors& ipol);
void calc_energy_shift_gamma_tabulated(extendedSourceData& dat_ind_a, const int ind_h, const int ind_x,
                                       double**** energy_shift_bins, const int num_bins, const int num_shift_bins,
                                       const double gamma);
double calc_energy_shift_gamma(extendedSourceTable* tab, const InterpFactors ipol,
                               const double bin_val, const double gamma);
static lpReflFrac* calc_refl_frac_ext(const extendedSourceTable* tab, const InterpFactors ipol, int* status);

double calc_primary_lensing(const double incl, const extendedSourceTable* prim_tab, const InterpFactors ipol);
void calc_ring_emis(emisProfile* emisProf, double a, double h, double x, double gamma, extendedSourceTable* tab, int* status);
void apply_flux_correction(extendedSourceData* dat, const double* r,
                           const double spin, const int n_h, const int n_x, const int n_rad);
void calc_ring_fractions(emisProfile* emisProf, const relParam* param, extendedSourceTable* tab, int* status);

// functions used globally:
void calc_emis_ring_source(emisProfile* emisProf, const relParam* param, int* status);
void calc_emis_slab_source(emisProfile* emisProf, const relParam* param, int* status);

double get_energyboost_ext_source_obs(const relParam* param, const double gamma);
double get_energyboost_ext_source_disk(const relParam* param, const double radius, const double gamma);

// constructors and destructors:
extendedSourceTable* new_extended_source_table(const int n_h, const int n_x,
                                               const int n_a, const int n_bin, const int n_storage,
                                               int* status);

extendedSourceData* new_extended_source_data(int* status, const bool is_lensing_table);

void free_extended_source_table(extendedSourceTable* tab, const bool is_lensing_table);

void free_extended_source_data(extendedSourceData* dat, const bool is_lensing_table);


#endif /* RELLP_EXTENDED_H_ */
