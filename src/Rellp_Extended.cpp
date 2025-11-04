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


/**
 * @brief This module creates emissivity profiles from extended sources ready for usage in relxill
 * @details This module reads values from extended tables produced with ray-tracing calculations
 * originally by Stefan Licklederer, later significantly reworked by Alexey Nekrasov
 */


#include "cfitsio.h"
#include "Rellp_Extended.h"
#include "Rellp.h"
#include "Relphysics.h"

extern "C"{
#include "relutility.h"
#include "writeOutfiles.h"
}

extendedSourceTable *cached_extendedSource_table = nullptr; // emissivity table

extendedSourceTable *cached_extendedSource_table_prim = nullptr; // lensing table

/**
 * @brief Gets extended source table with a given name EXT_TABLE_FILENAME
 * @details
 * @param status pointer, status of the code for debugging
 * @return cached_extendedSource_table
 */
extendedSourceTable* get_es_table(int* status){
  CHECK_STATUS_RET(*status, nullptr);

  if (cached_extendedSource_table == nullptr) {
    read_extendedSource_table(EXT_TABLE_FILENAME, &cached_extendedSource_table, status);
    CHECK_STATUS_RET(*status, nullptr);
  }
  return cached_extendedSource_table;
}

extendedSourceTable* get_prim_table(int* status){
  CHECK_STATUS_RET(*status, nullptr);

  bool is_lensing_table = true; // quick solution to have two tables instead of one
  if (cached_extendedSource_table_prim == nullptr) {
    read_extendedSource_table(PRIM_TABLE_FILENAME, &cached_extendedSource_table_prim, status,
                              is_lensing_table);
    CHECK_STATUS_RET(*status, nullptr);
  }
  return cached_extendedSource_table_prim;
}

/**
 * @brief Gets rows of the table corresponding to different spins
 * @details Uses cfitsio routines to read preloaded fits file
 * @param nrows integer, number of rows in the table
 * @param val two-fold pointer to the value which stores the table rows
 * @param fptr pointer to the file
 * @param status pointer status of the code for debugging
 */
static void get_spin_rows_extended_table(const int nrows, float **val, fitsfile *fptr, int *status){

  *val = (float*)  malloc(nrows* sizeof(float));

  for (int ii = 0; ii < nrows; ii++){

    int index = ii + 2; // number of the row we read in the fits table
    fits_movabs_hdu(fptr, index, nullptr, status);

    fits_read_key(fptr, TFLOAT, "SPIN", &((*val)[ii]), nullptr, status);
    if (*status != EXIT_SUCCESS) {
      CHECK_RELXILL_ERROR("reading spin dimension in table failed,", status);
    }
  }

}

/**
 * @brief Interpolates a given tabulated value (between 4 base points) between h and x,
 * and then between spins. In fact, this is 3d interpolation in (a,h,x)
 * @return interpolated value, double
 */
double get_interpolated_value_a_h_x(const double ifac_a, const double ifac_h[2], const double ifac_x[2],
                              const double vals_lo00, const double vals_lo01, const double vals_lo10, const double vals_lo11,
                                    const double vals_hi00, const double vals_hi01, const double vals_hi10, const double vals_hi11){
  // does h-r interpolation for spin[ind_a]
  double val_at_spin_lo = interp_lin_2d(ifac_x[0],ifac_h[0],
                                        vals_lo00,vals_lo01,vals_lo10,vals_lo11);
  // does h-r interpolation for spin[ind_a+1]
  double val_at_spin_hi = interp_lin_2d(ifac_x[1],ifac_h[1],
                                        vals_hi00,vals_hi01,vals_hi10,vals_hi11);
  // grid points must be equal, no adjustment done here
  return interp_lin_1d(ifac_a, val_at_spin_lo, val_at_spin_hi);
}

/**
 * @brief Divides flux by Kerr-to-flat area ratio and by Lorentz factor
 * TODO: cache these corrections in new table, so this function will be obsolete
 */
void apply_flux_correction(extendedSourceData *dat, double *r, double spin, int n_h, int n_x, int n_rad){
  for (int i=0; i<n_h; i++) {
    for (int j=0; j<n_x; j++) {
      for (int k=0; k<n_rad; k++){
        dat->val[i][j][k] /= (kerr_to_flat_area_ratio(r[k], spin) * gamma_correction(r[k], spin));
      }
    }
  }
}

/**
 * @brief Gets data for a single source from the table and then makes a relativistic correction of the flux
 * @details
 * @param fptr pointer, to the fits file location
 * @param n_h integer, number of heights simulated
 * @param n_x integer, number of ring radii simulated
 * @param n_bin integer, number of radial bins
 * @param n_storage integer, number of energyshift storage bins. Used only to read the fifth column with shift
 * @param spin double, spin value of the source
 * @return extendedSourceData, data from the table
 */
static extendedSourceData *
load_single_extendedSourceData(fitsfile *fptr, int n_h, int n_x,
                               int n_bin, int n_storage, int *status,
                               double spin, bool is_lensing_table) {

  extendedSourceData *dat = new_extendedSourceData(status, is_lensing_table);
  CHECK_MALLOC_RET_STATUS(dat, status, nullptr)

  int anynullptr = 0;
  double doublenullptr = 0.0;

  /*  read the columns */
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, n_bin, &doublenullptr,
                dat->bin,&anynullptr, status);

  for (int ii = 0; ii < n_h; ii++) {
    fits_read_col(fptr, TDOUBLE, 2, 1 + ii * n_x, 1, 1, &doublenullptr,
                  &(dat->height[ii]), &anynullptr, status);
    fits_read_col(fptr, TDOUBLE, 3, 1 + ii * n_x, 1, n_x,
                  &doublenullptr, dat->x[ii], &anynullptr, status);
  }

  for (int ii = 0; ii < n_h; ii++) {
    for (int jj = 0; jj < n_x; jj++) {

      int irow = (ii * n_x) + jj + 1;

      fits_read_col(fptr, TDOUBLE, 4, irow, 1, n_bin, &doublenullptr,
                    dat->val[ii][jj],&anynullptr, status);

      for (int kk = 0; kk < n_bin; kk++){
        // should never happen with the proper table:
        if (dat->val[ii][jj][kk] < 0.0) {
          CHECK_RELXILL_ERROR("Negative value loaded from table, check your table", status);
        }
        fits_read_col(fptr, TDOUBLE, 5, irow, kk*n_storage+1, n_storage,
                      &doublenullptr, dat->energy_data_storage[ii][jj][kk],&anynullptr, status);
      }

      // n_rad to 1 because it is not an array just a single value
      fits_read_col(fptr, TDOUBLE, 6, irow, 1, 1, &doublenullptr,
                    &dat->f_ad[ii][jj],&anynullptr, status);

      fits_read_col(fptr, TDOUBLE, 7, irow, 1, 1, &doublenullptr,
                    &dat->f_inf[ii][jj],&anynullptr, status);
      if (is_lensing_table) { // extra column which is now only in the primary table
        fits_read_col(fptr, TDOUBLE, 8, irow, 1, 1, &doublenullptr,
                      &dat->f_bh[ii][jj],&anynullptr, status);
      }
      dat->refl_frac[ii][jj] = dat->f_ad[ii][jj]/dat->f_inf[ii][jj];

    }
  }
  if (!is_lensing_table) { // only for flux, not for lensing_factor
    apply_flux_correction(dat, dat->bin, spin, n_h, n_x, n_bin);
  }
  return dat;
}


void get_ipol_factors_a_h_x(const double a, const double h, const double x,
                            const extendedSourceTable *tab, InterpFactors& ipol) {
  // get interpolation factors between spins
  get_ipol_factor((float) a, tab->a, tab->n_a, &ipol.ind_a, &ipol.ifac_a);
  if (ipol.ifac_a < 0.0 or ipol.ifac_a > 1.0){
  }
  for (int ii = 0; ii < 2; ii++) {
    // get interpolation factors between heights
    get_ipol_factor_double(h, tab->dat[ipol.ind_a + ii]->height, tab->n_h,
                           &(ipol.ind_h[ii]),&(ipol.ifac_h[ii]));
    if (ipol.ifac_h[ii] < 0.0 or ipol.ifac_h[ii] > 1.0){
      printf(" *** Warning : height interpolation becomes extrapolation: ifac_h = %.1f \n", ipol.ifac_h[ii]);
      printf(" *** The values: a = %.3f, h = %.3f, x = %.3f \n", a, h, x);
      printf("Better limit a range of values within table limits \n");
    }
    // interpolation factors between extents
    int index = (ipol.ifac_h[ii] > 0.999999) ? ipol.ind_h[ii] + 1 : ipol.ind_h[ii];
    // this additional index shift doesn't matter for the current table
    // but if tabulated extents will differ for different heights - it will matter
    get_ipol_factor_double(x, tab->dat[ipol.ind_a + ii]->x[index], tab->n_x,
                           &(ipol.ind_x[ii]), &(ipol.ifac_x[ii]));
    if (ipol.ifac_x[ii] < 0.0 or ipol.ifac_x[ii] > 1.0){
      printf(" *** Warning : radius interpolation becomes extrapolation: ifac_x = %.1f \n", ipol.ifac_x[ii]);
      printf("Better limit a range of values within table limits \n");
    }
  }
}


/**
 * @brief determine the gshift^Gamma for table grid values. Required for further interpolation
 */
void calc_energy_shift_gamma_tabulated(extendedSourceData &dat_ind_a, const int ind_h, const int ind_x,
                                       double ****energy_shift_bins, int num_bins, int num_shift_bins, double gamma){
  for (int i_h : {ind_h, ind_h + 1}) { // vertical bins, ind_h
    for (int j_x : {ind_x, ind_x + 1}) { // radial bins, ind_x
      for (int k_b = 0; k_b < num_bins; k_b++) { // disk radial bins or lensing inclination bins
        double _sum = 0.0; // sum is zeroed every time
        for (int l_e = 0; l_e < num_shift_bins; l_e++) { // energy shift bins
          _sum += pow(energy_shift_bins[i_h][j_x][k_b][l_e], gamma);
          assert(_sum >= 0.0 || !"summation of g is negative, wrong value in table\n");
        }
        dat_ind_a.gshift[i_h][j_x][k_b] = _sum / num_shift_bins;
      }
    }
  }
}

/**
 * @details The difference to calc_energy_shift_gamma_tabulated function is that here we are not on the grid points
 * of the table, so we need to do extra actions. This function is more general than calc_energy_shift_gamma_tabulated
 */
double calc_energy_shift_gamma(const extendedSourceTable *tab, InterpFactors ipol, double bin_val, double gamma) {
  int ind_bin_lo, ind_bin_hi;
  double ifac_bin_lo, ifac_bin_hi;
  get_ipol_factor_double(bin_val, tab->dat[ipol.ind_a]->bin, tab->n_bin, &ind_bin_lo, &ifac_bin_lo);
  get_ipol_factor_double(bin_val, tab->dat[ipol.ind_a + 1]->bin, tab->n_bin, &ind_bin_hi, &ifac_bin_hi);

  for (int i_h : {ipol.ind_h[0], ipol.ind_h[0] + 1}) { // vertical bins
    for (int j_x : {ipol.ind_x[0], ipol.ind_x[0] + 1}) { // radial extent bins
      for (int k_b : {ind_bin_lo, ind_bin_lo + 1}) { // disk radial bins or lensing inclination bins
        double _sum = 0.0; // sum is zeroed every time
        for (int l_e = 0; l_e < tab->n_storage; l_e++) { // energy shift bins
          _sum += pow(tab->dat[ipol.ind_a]->energy_data_storage[i_h][j_x][k_b][l_e], gamma);
          assert(_sum >= 0.0 || !"summation of g is negative, wrong value in table\n");
        }
        tab->dat[ipol.ind_a]->gshift[i_h][j_x][k_b] = _sum / tab->n_storage;
      }
    }
  }
  for (int i_h : {ipol.ind_h[1], ipol.ind_h[1] + 1}) { // vertical bins
    for (int j_x : {ipol.ind_x[1], ipol.ind_x[1] + 1}) { // radial extent bins
      for (int k_b : {ind_bin_hi, ind_bin_hi + 1}) {
        double _sum = 0.0;
        for (int l_e = 0; l_e < tab->n_storage; l_e++) {
          _sum += pow(tab->dat[ipol.ind_a + 1]->energy_data_storage[i_h][j_x][k_b][l_e], gamma);
          assert(_sum >= 0.0 || !"summation of g is negative, wrong value in table\n");
        }
        tab->dat[ipol.ind_a + 1]->gshift[i_h][j_x][k_b] = _sum / tab->n_storage;
      }
    }
  }

  double tab_lo_00 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_bin_lo],
                                   tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_bin_lo + 1]);
  double tab_lo_01 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_bin_lo],
                                   tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_bin_lo + 1]);
  double tab_lo_10 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_bin_lo],
                                   tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_bin_lo + 1]);
  double tab_lo_11 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_bin_lo],
                                   tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_bin_lo + 1]);

  double tab_hi_00 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_bin_hi],
                                   tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_bin_hi + 1]);
  double tab_hi_01 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_bin_hi],
                                   tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_bin_hi + 1]);
  double tab_hi_10 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_bin_hi],
                                   tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_bin_hi + 1]);
  double tab_hi_11 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_bin_hi],
                                   tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_bin_hi + 1]);

  double g_gamma = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                tab_lo_00, tab_lo_01, tab_lo_10, tab_lo_11,
                                                tab_hi_00, tab_hi_01, tab_hi_10, tab_hi_11);
  assert (g_gamma >= 0.0 || !"g^Gamma negative, wrong interpolation!\n");
  return g_gamma;
}

/**
 * @brief Interpolate the table between two spins
 * @details
 * @param tab; extendedSourceTable, pointer to the table as declared above
 * @param ipol; structure with interpolation factors between a, h, x
 * @param gamma; double, photon index
 * @return emis_profile_table; emisProfile, pointer to the structure with emissivity profile
 */
emisProfile *interpol_extended_table(extendedSourceTable *tab, InterpFactors ipol, double gamma, double* radial_grid, int *status) {
  CHECK_STATUS_RET(*status, nullptr);

  // The radial grid of the table is different for different spins (risco(a)), so we must interpolate it
  assert(radial_grid == nullptr);
  radial_grid = (double*) malloc(sizeof(double)*tab->n_bin);
  CHECK_MALLOC_RET_STATUS(radial_grid,status,nullptr)
  for (int ii = 0; ii < tab->n_bin; ii++) {
    radial_grid[ii] = interp_lin_1d(ipol.ifac_a, tab->dat[ipol.ind_a]->bin[ii], tab->dat[ipol.ind_a+1]->bin[ii]);
  }

  emisProfile* emis_profile_table = new_emisProfile(radial_grid, tab->n_bin, status);

  calc_energy_shift_gamma_tabulated(*tab->dat[ipol.ind_a], ipol.ind_h[0], ipol.ind_x[0],
                                    tab->dat[ipol.ind_a]->energy_data_storage,
                          tab->n_bin, tab->n_storage, gamma);
  calc_energy_shift_gamma_tabulated(*tab->dat[ipol.ind_a + 1], ipol.ind_h[1], ipol.ind_x[1],
                                    tab->dat[ipol.ind_a + 1]->energy_data_storage,
                          tab->n_bin, tab->n_storage, gamma);

  // pre-allocation of values shifted in radial grid before interpolation
  double tab_lo_00, tab_lo_01, tab_lo_10, tab_lo_11, tab_hi_00, tab_hi_01, tab_hi_10, tab_hi_11;
  double fac_rlo, fac_rhi;
  int ind_rlo, ind_rhi;
  for (int ii = 0; ii < tab->n_bin; ii++) {
    // As radial grids shift between spins, we try to compensate for this shift in the values
    get_ipol_factor_double(radial_grid[ii], tab->dat[ipol.ind_a]->bin, tab->n_bin, &ind_rlo, &fac_rlo);
    get_ipol_factor_double(radial_grid[ii], tab->dat[ipol.ind_a + 1]->bin, tab->n_bin, &ind_rhi, &fac_rhi);

    // Interpolation of the flux
    tab_lo_00 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo],
                              tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_01 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo],
                              tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo + 1]);
    tab_lo_10 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo],
                              tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_11 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo],
                              tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo + 1]);

    tab_hi_00 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_01 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi + 1]);
    tab_hi_10 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_11 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi + 1]);

    emis_profile_table->emis[ii] = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                                tab_lo_00, tab_lo_01,
                                                                tab_lo_10, tab_lo_11,
                                                                tab_hi_00, tab_hi_01,
                                                                tab_hi_10, tab_hi_11);

    // Re-use temporary variables for both emis and g_gamma
    tab_lo_00 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo],
                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_01 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo],
                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo + 1]);
    tab_lo_10 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo],
                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_11 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo],
                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo + 1]);

    tab_hi_00 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_01 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi + 1]);
    tab_hi_10 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_11 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi],
                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi + 1]);

    double g_gamma = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                  tab_lo_00, tab_lo_01,
                                                  tab_lo_10, tab_lo_11,
                                                  tab_hi_00, tab_hi_01,
                                                  tab_hi_10, tab_hi_11);

    emis_profile_table->emis[ii] *= g_gamma;
    if (emis_profile_table->emis[ii] < 0.0) {
      CHECK_RELXILL_ERROR("Negative emissivity bin caught after table interpolation!", status);
    }
  }
  return emis_profile_table;
}


/**
 * @brief Calculates reflection fraction from interpolated table values, and photon fractions
 * @details Can it be combined/restructured with interpol_extended_table ?
 * @param param, const relParam, a pointer to the structure of parameters used, the parameters include:
 * a - spin, h - height, x - ring radius
 * @param tab, extendedSourceTable, pointer to the structure containing a table
 * @return lpReflFrac, pointer to the structure containing interpolated reflection fraction and photon fractions
 */
static lpReflFrac *calc_refl_frac_ext(const extendedSourceTable *tab, InterpFactors ipol, int *status) {

  lpReflFrac *fractions = new_lpReflFrac(status);
  CHECK_STATUS_RET(*status, fractions);

  fractions->f_ad = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                 tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0]][ipol.ind_x[0]],
                                                 tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0]][ipol.ind_x[0] + 1],
                                                 tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0] + 1][ipol.ind_x[0]],
                                                 tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1],
                                                 tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1]][ipol.ind_x[1]],
                                                 tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1]][ipol.ind_x[1] + 1],
                                                 tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1] + 1][ipol.ind_x[1]],
                                                 tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1]);

  fractions->f_inf = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                  tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0]][ipol.ind_x[0]],
                                                  tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0]][ipol.ind_x[0] + 1],
                                                  tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0] + 1][ipol.ind_x[0]],
                                                  tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1],
                                                  tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1]][ipol.ind_x[1]],
                                                  tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1]][ipol.ind_x[1] + 1],
                                                  tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1] + 1][ipol.ind_x[1]],
                                                  tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1]);

  fractions->f_bh = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                 tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0]][ipol.ind_x[0]],
                                                 tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0]][ipol.ind_x[0] + 1],
                                                 tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0] + 1][ipol.ind_x[0]],
                                                 tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1],
                                                 tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1]][ipol.ind_x[1]],
                                                 tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1]][ipol.ind_x[1] + 1],
                                                 tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1] + 1][ipol.ind_x[1]],
                                                 tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1]);

  fractions->refl_frac = fractions->f_ad / fractions->f_inf;

  if (fractions->f_inf > 0.5) { // photons are not allowed to cross the disk plane - should never be the case
    fractions->f_inf = 0.5;
    printf("Rellp_Extended: escaping photon fraction more than 0.5 (happens only for h>40rg). \n");
    printf("Replacing with 0.5. This will be fixed in the new extended table. \n");
  }

  // these are equal for now. f_inf_rest is used to calculate primary spectrum
  fractions->f_inf_rest = fractions->f_inf; // true as long as beta=0

  return fractions;
}

/*
 * Main routines to be used in the project
 * */
/**
 * @brief Gets ring emissivities and normalization factors to emisProfile structure
 * @details Calculates emissivity for the extended source using calc_ring_emis function
 * @param emisProf, pointer to the structure containing emissivity profile
 * @param param, pointer to the structure containing the parameters used
 */
void calc_emis_ring_source(emisProfile *emisProf, const relParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  extendedSourceTable* emissivity_table = get_es_table(status);
  extendedSourceTable* prim_table = get_prim_table(status);

  calc_ring_emis(emisProf, param->a, param->height, param->x, param->gamma, emissivity_table, status);
  calc_ring_fractions(emisProf, param, prim_table, status);

  CHECK_STATUS_VOID(*status);

}

void calc_emis_slab_source(emisProfile *emisProf, const relParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  extendedSourceTable* emissivity_table = get_es_table(status);
  extendedSourceTable* prim_table = get_prim_table(status);
  InterpFactors ipol;
  InterpFactors prim_ipol;
  // TODO more advanced scheme for primary normalization, currently assume just the values from the outermost ring:
  calc_ring_fractions(emisProf, param, prim_table, status);

  const int num_rings = 100;
  std::array<double, num_rings> slab_grid{};
  emisProfile *slab_emis_sum = nullptr;
  double* rad_grid = nullptr;

  // auto *param_tmp = const_cast<relParam *>(param); // it could work too
  double weight_i;
  double weight_sum = 0.0;

  for (int ii = 0; ii < num_rings; ii++) { // define grid of rings for slab
    slab_grid[ii] = param->x_in + (param->x - param->x_in) * ii / (num_rings - 1);
  }
  for (int ii = 0; ii < num_rings; ii++) {

    get_ipol_factors_a_h_x(param->a, param->height, slab_grid[ii], emissivity_table, ipol);
    get_ipol_factors_a_h_x(param->a, param->height, slab_grid[ii], prim_table, prim_ipol);

    emisProfile *ring_emis_tmp = interpol_extended_table(emissivity_table, ipol, param->gamma, rad_grid, status);
    if (ii == 0) {
      slab_emis_sum = new_emisProfile(ring_emis_tmp->re, ring_emis_tmp->nr, status);
    }
    assert(slab_emis_sum != nullptr);
    assert(slab_emis_sum->emis != nullptr);

    ring_emis_tmp->photon_fate_fractions = calc_refl_frac_ext(prim_table, prim_ipol, status);

    if (ii == num_rings - 1) { // last bin
      weight_i = slab_grid[num_rings - 1] * (slab_grid[num_rings - 1] -
          slab_grid[num_rings - 2]) * ring_emis_tmp->photon_fate_fractions->refl_frac;
    } else { // ordinary bin
      weight_i = slab_grid[ii] * (slab_grid[ii + 1
                                            ] -
          slab_grid[ii]) * ring_emis_tmp->photon_fate_fractions->refl_frac;
    }
    weight_sum += weight_i;
    for (int nn = 0; nn < ring_emis_tmp->nr; nn++) {
      slab_emis_sum->emis[nn] += ring_emis_tmp->emis[nn] * weight_i;
    }
    // free(ring_emis_tmp->re);
    // free_emisProfile(ring_emis_tmp);  // TODO: why not?? I think you definitely need to free it [TD]
  }
  for (int nn = 0; nn < slab_emis_sum->nr; nn++) {
    slab_emis_sum->emis[nn] /= weight_sum; // normalize final profile with weights (one could also divide by disk area)
  }
  rebin_emisprofile_on_radial_grid(emisProf, slab_emis_sum, status);

  free(rad_grid);
  free_emisProfile(slab_emis_sum);

  CHECK_STATUS_VOID(*status);
}

/**
 * @brief Secondary main routine which gets g^Gamma from source to observer from the table (at relxill inclination)
 * @details for now it is almost duplicated with the function below, which does the same for the disk point
 * @return The g^Gamma value with Gamma not necessarily taken from relParam (e.g. can be 1 for just energyshift)
 */
double get_energyboost_ext_source_obs(const relParam *param, double gamma) {
  int status = EXIT_SUCCESS; // in case it is loaded from somewhere in the code
  extendedSourceTable* prim_table = get_prim_table(&status);
  InterpFactors prim_ipol;
  get_ipol_factors_a_h_x(param->a, param->height, param->x, prim_table, prim_ipol);
  // input incl can be both in radian and in grads. The table has cos incl, so convert it here to cos
  // (I use the fact that relxill limits are 3 and 87 deg, so if input incl is < 1.57 (~pi/2) - it is radian)
  double cos_incl = (param->incl < 0.5 * M_PI) ? cos(param->incl) : cos(param->incl * M_PI / 180.0);
  double gshift = calc_energy_shift_gamma(prim_table, prim_ipol, cos_incl, gamma);
  return gshift;
}

/**
 * @brief Secondary main routine which gets g^Gamma from source to disk for a given disk radius
 * @param param relxill parameters and values
 * @param radius Disk radius
 * @param gamma Optionally, gamma can be different from the gamma spectrum index (e.g., 1, when we need just energy shift)
 * @return The g^Gamma value
 */
double get_energyboost_ext_source_disk(const relParam *param, double radius, double gamma) {
  int status = EXIT_SUCCESS;
  extendedSourceTable* table = get_es_table(&status);
  InterpFactors ipol;
  get_ipol_factors_a_h_x(param->a, param->height, param->x, table, ipol);
  double gshift = calc_energy_shift_gamma(table, ipol, radius, gamma);
  return gshift;
}


/**
 * @brief Reads extended source table
 * @details
 * @param filename; const char, pointer to the name of the file (declared in Rellp_Extended.h)
 * @param inp_tab; extendedSourceTable, two-fold pointer to the input table data
 */
void read_extendedSource_table(const char *filename, extendedSourceTable **inp_tab, int *status, bool is_lensing_table) {

  extendedSourceTable *tab = (*inp_tab);
  fitsfile *fptr = nullptr;

  char fullfilename[999];

  do {
    if (tab != nullptr) {
      RELXILL_ERROR("extended source table already loaded\n", status);
      break;
    }
    if (is_lensing_table) {
      tab = new_extendedSourceTable(PRIM_TABLE_NHEIGHT, PRIM_TABLE_NDISTANCE, PRIM_TABLE_NA,
                                    PRIM_TABLE_NINCL, PRIM_TABLE_NENERGYBINS, status);
    } else {
      tab = new_extendedSourceTable(EXT_TABLE_NHEIGHT, EXT_TABLE_NDISTANCE, EXT_TABLE_NA,
                                    EXT_TABLE_NRAD, EXT_TABLE_NENERGYBINS, status);
    }

    CHECK_STATUS_BREAK(*status);

    // get the full filename
    if (sprintf(fullfilename, "%s/%s", get_relxill_table_path(), filename) == -1) {
      RELXILL_ERROR("failed to construct full path the rel table\n", status);


      break;
    }

    // open the file
    if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
      CHECK_RELXILL_ERROR("opening of the extended source table failed", status);
      printf("    full path given: %s \n", fullfilename);
      break;
    }

    get_spin_rows_extended_table(tab->n_a, &(tab->a), fptr, status);
    CHECK_RELXILL_ERROR("reading of spin axis failed\n", status);

    // now load the full table (need to go through all extensions)
    int ii;
    int rownum;
    for (ii = 0; ii < tab->n_a; ii++) {

      rownum = ii + 2; // number of the row we read in the fits table
      fits_movabs_hdu(fptr, rownum, nullptr, status);

      if (tab->dat[ii] != nullptr) {
        CHECK_RELXILL_ERROR("reading of spin axis failed\n", status);
      }

      tab->dat[ii] = load_single_extendedSourceData(
          fptr, tab->n_h, tab->n_x, tab->n_bin, tab->n_storage, status,
          tab->a[ii], is_lensing_table);
      if (*status != EXIT_SUCCESS) {
        RELXILL_ERROR("failed to load data from the extended source table into memory\n", status);

        break;
      }
    }

  } while (false);

  if (*status == EXIT_SUCCESS) {
    (*inp_tab) = tab;
  } else {
    free_extendedSourceTable(tab, is_lensing_table);
  }

  if(fptr != nullptr) {fits_close_file(fptr, status);}

}

double calc_primary_lensing(double incl, extendedSourceTable *prim_tab, InterpFactors ipol) {
  int ind_incl;
  double ifac_incl;
  get_ipol_factor_double(incl, prim_tab->dat[ipol.ind_a]->bin, prim_tab->n_bin,
                         &ind_incl, &ifac_incl);
  if (ifac_incl < 0.0 or ifac_incl > 1.0){
    // should not happen in the current table
    printf(" *** Warning : inclination interpolation becomes extrapolation: ifac_i = %.1f \n", ifac_incl);
    printf("Better limit a range of values within table limits, 0.05 < cos(incl) < 1.0 \n");
  }

  double lens_low = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                 prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][ind_incl],
                                                 prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] +
                                                                                               1][ind_incl],
                                                 prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] +
                                                                                1][ipol.ind_x[0]][ind_incl],
                                                 prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] +
                                                                                                   1][ind_incl],
                                                 prim_tab->dat[ipol.ind_a +
                                                               1]->val[ipol.ind_h[1]][ipol.ind_x[1]][ind_incl],
                                                 prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] +
                                                                                                   1][ind_incl],
                                                 prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] +
                                                                                    1][ipol.ind_x[1]][ind_incl],
                                                 prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] +
                                                                                                       1][ind_incl]);
  double lens_high = get_interpolated_value_a_h_x(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                  prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][
                                                      ind_incl + 1],
                                                  prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] + 1][
                                                      ind_incl + 1],
                                                  prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0]][
                                                      ind_incl + 1],
                                                  prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][
                                                      ind_incl + 1],
                                                  prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1]][
                                                      ind_incl + 1],
                                                  prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] + 1][
                                                      ind_incl + 1],
                                                  prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1]][
                                                      ind_incl + 1],
                                                  prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] +
                                                                                                        1][ind_incl +
                                                                                                           1]);
  double lensing = interp_lin_1d(ifac_incl, lens_low, lens_high);
  return lensing;
}


/**
 * @brief Calculates ring emissivity
 * @details
 * @param emisProf; emisProfile, pointer to the structure with emissivity profile
 * @param param; relParam, pointer to the parameters structure
 * @param tab; extendedSourceTable, pointer to the table structure
 */
void calc_ring_emis(emisProfile *emisProf, const double a, const double h, const double x, const double gamma,
                    extendedSourceTable *tab, int *status) {

  CHECK_STATUS_VOID(*status);

  InterpFactors ipol;
  get_ipol_factors_a_h_x(a, h, x, tab, ipol);

  double* radial_grid = nullptr;
  emisProfile *emis_profile_table = interpol_extended_table(tab, ipol, gamma, radial_grid, status);
  rebin_emisprofile_on_radial_grid(emisProf, emis_profile_table, status);

  free(radial_grid); // radial_grid is _not_ owned by emissivity profile structure
  free_emisProfile(emis_profile_table);


  // TODO take the beaming of the jet into account (see Dauser et al., 2013)
//  for (int ii=0; ii<tab->n_rad; ii++){
//    if (param->beta > 1e-6) {
      // emisProf->emis[ii] *= pow(doppler_factor(emisProf->del_emit[ii], param->beta), 2);
//    }
//  }
}

/**
 * @brief Calculates photon fractions (esc, bh, disk), lensing and boost of the primary emission for ring geometry
 * @details
 * @param emisProf; emisProfile, pointer to the structure with emissivity profile and necessary quantities
 * @param param; relParam, pointer to the parameters structure
 * @param prim_tab; extendedSourceTable, pointer to the primary table
 */
void calc_ring_fractions(emisProfile *emisProf, const relParam *param, extendedSourceTable *prim_tab, int *status) {
  CHECK_STATUS_VOID(*status);

  InterpFactors prim_ipol;
  get_ipol_factors_a_h_x(param->a, param->height, param->x, prim_tab, prim_ipol);
  emisProf->photon_fate_fractions = calc_refl_frac_ext(prim_tab, prim_ipol, status);
  // input incl can be both in radian and in grads. The table has cos incl, so convert:
  double cos_incl = (param->incl < 0.5 * M_PI) ? cos(param->incl) : cos(param->incl * M_PI / 180.0);
  emisProf->photon_fate_fractions->lensing = calc_primary_lensing(cos_incl, prim_tab, prim_ipol);
  emisProf->photon_fate_fractions->lensing_and_boost_factor = (emisProf->photon_fate_fractions->lensing *
      calc_energy_shift_gamma(prim_tab, prim_ipol, cos_incl, param->gamma));
}


// constructors and destructors //
/**
 *
 * @brief Sets the memory for new extended source table
 * @details
 * @param n_h; int, number of heights in the table
 * @param n_r; int, number of ring sources in the table
 * @param n_a; int, number of spins in the table
 * @param n_bin; int, number of bins in the table (radius or inclination bins)
 * @param n_storage; int, number of the energy shift storage in the table
 * @return new_extendedSourceTable, new extended source table
 */
extendedSourceTable *new_extendedSourceTable(int n_h, int n_r,
                                             int n_a, int n_bin, int n_storage,
                                             int *status) {
  auto *tab = (extendedSourceTable *) malloc(sizeof(extendedSourceTable));

  CHECK_MALLOC_RET_STATUS(tab, status, nullptr)

  tab->n_a = n_a;
  tab->n_h = n_h;
  tab->n_bin = n_bin;
  tab->n_storage = n_storage;
  tab->n_x = n_r;

  tab->a = nullptr;

  tab->dat = nullptr;

  tab->dat = (extendedSourceData **) malloc(sizeof(extendedSourceData *) * tab->n_a);

  CHECK_MALLOC_RET_STATUS(tab->dat, status, tab)

  int ii;
  for (ii = 0; ii < tab->n_a; ii++) {
    tab->dat[ii] = nullptr;
  }
  return tab;
}


/**
 *
 * @brief Sets the memory for new extended source data
 * @details
 * @return extendedSourceData, new extended source data
 */
extendedSourceData *new_extendedSourceData(int *status, bool is_lensing_table) {
  int n_h = EXT_TABLE_NHEIGHT;
  int n_x = EXT_TABLE_NDISTANCE;
  int n_bin = EXT_TABLE_NRAD;
  int n_shift = EXT_TABLE_NENERGYBINS;
  if (is_lensing_table) {
    n_h = PRIM_TABLE_NHEIGHT;
    n_x = PRIM_TABLE_NDISTANCE;
    n_bin = PRIM_TABLE_NINCL;
    n_shift = PRIM_TABLE_NENERGYBINS;
  }
  auto *dat = (extendedSourceData *) malloc(sizeof(extendedSourceData));
  CHECK_MALLOC_RET_STATUS(dat, status, nullptr)

  dat->height = (double *) malloc(sizeof(double) * n_h);

  CHECK_MALLOC_RET_STATUS(dat->height, status, nullptr)

  dat->x = (double **) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->x, status, nullptr)

  dat->bin = (double *) malloc(sizeof(double) * n_bin);
  CHECK_MALLOC_RET_STATUS(dat->bin, status, nullptr)

  dat->refl_frac = (double **) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, nullptr)

  dat->f_ad = (double **) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->f_ad, status, nullptr)

  dat->f_inf = (double **) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->f_inf, status, nullptr)

  dat->f_bh = (double **) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->f_bh, status, nullptr)

  dat->val = (double ***) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->val, status, nullptr)

  dat->gshift = (double ***) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->gshift, status, nullptr)

  dat->energy_data_storage = (double ****) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, nullptr)

  for (int i=0; i < n_h; i++){

    dat->refl_frac[i] = (double *) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, nullptr)

    dat->x[i] = (double *) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->x, status, nullptr)

    dat->f_ad[i] = (double *) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->f_ad, status, nullptr)

    dat->f_inf[i] = (double *) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->f_inf, status, nullptr)

    dat->f_bh[i] = (double *) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->f_bh, status, nullptr)

    dat->gshift[i] = (double **) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->gshift, status, nullptr)

    dat->val[i] = (double **) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->val, status, nullptr)

    dat->energy_data_storage[i] = (double ***) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, nullptr)

    for(int j=0; j<n_x; j++){
      dat->gshift[i][j] = (double *) malloc(sizeof(double) * n_bin);
      dat->val[i][j] = (double *) malloc(sizeof(double) * n_bin);
      dat->energy_data_storage[i][j] = (double **) malloc(sizeof(double) * n_bin);

      for(int k=0; k<n_bin; k++){
        dat->energy_data_storage[i][j][k] = (double *) malloc(sizeof(double) * n_shift);
      }
    }
  }

  return dat;
}


/**
 * @brief Frees the memory stored before for the data
 * @details
 * @param tab; extendedSourceTable, table to be freed from memory
 */
void free_extendedSourceTable(extendedSourceTable *tab, bool lensing) {
  if (tab != nullptr) {
    if (tab->dat != nullptr) {
      for (int ii = 0; ii < tab->n_a; ii++) {
        if (tab->dat[ii] != nullptr){
          free_extendedSourceData(tab->dat[ii], lensing);
          free(tab->dat[ii]);
        }
      }
      free(tab->dat);
    }
    free(tab->a);
    free(tab);
  }
}


/**
 * @brief Frees the memory stored before for the data
 * @details
 * @param dat; extendedSourceData, data to be freed
 */
void free_extendedSourceData(extendedSourceData *dat, bool is_lensing_table) {
  int n_h = EXT_TABLE_NHEIGHT;
  int n_x = EXT_TABLE_NDISTANCE;
  int n_r = EXT_TABLE_NRAD;
  if (is_lensing_table) {
    n_h = PRIM_TABLE_NHEIGHT; // re-define in case of lensing table, dimensions are different
    n_x = PRIM_TABLE_NDISTANCE;
    n_r = PRIM_TABLE_NINCL;
  }
  if (dat != nullptr) {
    for (int ii=0; ii < n_h; ii++) {
      if (dat->x != nullptr) free(dat->x[ii]);
      if (dat->val != nullptr) free(dat->val[ii]);
      if (dat->gshift != nullptr) free(dat->gshift[ii]);
      if (dat->refl_frac != nullptr) free(dat->refl_frac[ii]);
      if (dat->f_ad != nullptr) free(dat->f_ad[ii]);
      if (dat->f_inf != nullptr) free(dat->f_inf[ii]);
      if (dat->f_bh != nullptr) free(dat->f_bh[ii]);
      if (dat->energy_data_storage != nullptr) free(dat->energy_data_storage[ii]);

      for (int jj = 0; jj < n_x; jj++) {
        if (dat->val[ii] != nullptr) free(dat->val[ii][jj]);
        if (dat->gshift[ii] != nullptr) free(dat->gshift[ii][jj]);

        for (int kk = 0; kk < n_r; kk++){
          if (dat->energy_data_storage[ii][jj] != nullptr) free(dat->energy_data_storage[ii][jj][kk]);
        }
      }
    }
  }
    free(dat->bin);
    free(dat->height);
    free(dat->x);
    free(dat->val);
    free(dat->gshift);
    free(dat->energy_data_storage);
    free(dat->refl_frac);
    free(dat->f_inf);
    free(dat->f_ad);
    free(dat->f_bh);

}
