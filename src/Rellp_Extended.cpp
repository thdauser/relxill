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
 * @brief This module creates emissivity profile from the radially extended source ready for usage in relxill
 * @details This module reads values from the extended table produced with ray-tracing calculations
 * originally "corona-lean" by Stefan Licklederer, now significantly reworked by Alexey Nekrasov
 */


#include "cfitsio.h"
// #include "stdio.h"
//#include "Relbase.h"
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
 * @brief Gets extended source table with a given name
 * @details
 * @param status pointer, status of the code for debugging
 * @return extendedSourceTable, cached extended source table
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
  // maybe not the best decision, but otherwise how do I do all if else statements for two similar but different tables
  // without duplicating functions
  bool is_lensing_table = true;
  if (cached_extendedSource_table_prim == nullptr) {
    read_extendedSource_table(PRIM_TABLE_FILENAME, &cached_extendedSource_table_prim, status,
                              is_lensing_table);
    CHECK_STATUS_RET(*status, nullptr);
  }
  return cached_extendedSource_table_prim;
}

/**
 * @brief Gets spin values from the table for each row
 * @details Uses cfitsio to read preloaded fits file
 * @param nrows integer, number of rows in the table
 * @param val two-fold pointer, to the value to be read (better rename)
 * @param fptr pointer, to the file
 * @param status pointer, status of the code for debugging
 */
static void get_extendedsourcetable_spin(int nrows, float **val, fitsfile *fptr, int *status){

  *val = (float*)  malloc(nrows* sizeof(float));

  for (int i=0; i<nrows; i++){

    int index = i + 2; // number of the row we read in the fits table
    fits_movabs_hdu(fptr, index, nullptr, status);

    // dereference val to write a value to array
    fits_read_key(fptr, TFLOAT, "SPIN", &((*val)[i]), nullptr, status);
    assert(*status == EXIT_SUCCESS);
  }

}

/**
 * @brief This function interpolates a given tabulated value (between 4 base points) between h and r,
 * and then between spins. In fact, this is 3d interpolation in (a,h,x)
 * @details Probably the input can be a reference to the structure with values, and interp factors can be also struct
 * (Will it be safer, faster?)
 * @return The interpolated value, double
 */
double get_interpolated_value(const double ifac_a, const double ifac_h[2], const double ifac_r[2],
                              double vals_lo00, double vals_lo01, double vals_lo10, double vals_lo11,
                              double vals_hi00, double vals_hi01, double vals_hi10, double vals_hi11){
  // does h-r interpolation for spin[ind_a]
  double val_at_spin_lo = interp_lin_2d(ifac_r[0],ifac_h[0],
                                        vals_lo00,vals_lo01,vals_lo10,vals_lo11);
  // does h-r interpolation for spin[ind_a+1]
  double val_at_spin_hi = interp_lin_2d(ifac_r[1],ifac_h[1],
                                        vals_hi00,vals_hi01,vals_hi10,vals_hi11);
  // grid points must be equal, no adjustment done here
  return interp_lin_1d(ifac_a, val_at_spin_lo, val_at_spin_hi);
}

/**
 * @brief Simply divides flux by Kerr area (or kerr_flat ratio) and Lorentz correction
 * When I do this in table - this function is gone TODO: cache these corrections in new table
 */
void apply_flux_correction(extendedSourceData *dat, double *r, double spin, int n_h, int n_r, int n_rad){
  for (int i=0; i<n_h; i++) {
    for (int j=0; j<n_r; j++) {
      for (int k=0; k<n_rad; k++){
        dat->val[i][j][k] /= (kerr_to_flat_area_ratio(r[k], spin) * gamma_correction(r[k], spin));
      }
    }
  }
}

/**
 * @brief Gets data for a single source from the table and then makes a relativistic correction of the flux
 * @details Need to set the variables for fits_read_col instead of using numbers with unclear origin
 * Should the relativistic correction be made by a separate function better?
 * @param fptr pointer, to the fits file location
 * @param n_h integer, number of heights simulated
 * @param n_r integer, number of ring radii simulated
 * @param n_bin integer, number of radial bins
 * @param n_storage integer, number of energyshift storage bins. Used only to read the fifth column with shift
 * @param spin double, spin value of the source
 * @return extendedSourceData, data from the table
 */
static extendedSourceData *
load_single_extendedSourceData(fitsfile *fptr, int n_h, int n_r,
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
    fits_read_col(fptr, TDOUBLE, 2, 1 + ii * n_h, 1, 1, &doublenullptr,
                  &(dat->height[ii]),&anynullptr, status);
    fits_read_col(fptr, TDOUBLE, 3, 1 + ii * n_h, 1, n_r,
                  &doublenullptr, dat->r_ring[ii],&anynullptr, status);
  }

  for (int ii = 0; ii < n_h; ii++) {
    for (int jj = 0; jj < n_r; jj++) {

      int irow = (ii * n_r) + jj + 1;

      fits_read_col(fptr, TDOUBLE, 4, irow, 1, n_bin, &doublenullptr,
                    dat->val[ii][jj],&anynullptr, status);

      for (int kk = 0; kk < n_bin; kk++){
        assert(dat->val[ii][jj][kk] >= 0.0);
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
    /// TODO to be removed with new table, which does these corrections by default
    apply_flux_correction(dat, dat->bin, spin, n_h, n_r, n_bin);
  }
  return dat;
}


void get_ipol_factors_a_h_x(const relParam *param, const extendedSourceTable *tab, InterpFactors& ipol) {
  // interpolation factors between spins
  get_ipol_factor((float) param->a, tab->a, tab->n_a, &ipol.ind_a, &ipol.ifac_a);
  if (ipol.ifac_a < 0.0 or ipol.ifac_a > 1.0){
    // should not happen in the current table
    printf(" *** Warning : spin interpolation becomes extrapolation: ifac_a = %.1f \n", ipol.ifac_a);
    printf("Better limit a range of values within table limits \n");
  }
  for (int ii = 0; ii < 2; ii++) {
    // interpolation factors between heights
    get_ipol_factor_double(param->height, tab->dat[ipol.ind_a + ii]->height, tab->n_h,
                           &(ipol.ind_h[ii]),&(ipol.ifac_h[ii]));
    if (ipol.ifac_h[ii] < 0.0 or ipol.ifac_h[ii] > 1.0){
      // should not happen in the current table
      printf(" *** Warning : height interpolation becomes extrapolation: ifac_h = %.1f \n", ipol.ifac_h[ii]);
      printf("Better limit a range of values within table limits \n");
    }
    // interpolation factors between extents
    int index = (ipol.ifac_h[ii] > 0.999999) ? ipol.ind_h[ii] + 1 : ipol.ind_h[ii];
    // this additional index shift doesn't matter for the current table
    // but if tabulated extents will differ for different heights - it will matter
    get_ipol_factor_double(param->r_ring, tab->dat[ipol.ind_a + ii]->r_ring[index], tab->n_r,
                           &(ipol.ind_x[ii]), &(ipol.ifac_x[ii]));
    if (ipol.ifac_x[ii] < 0.0 or ipol.ifac_x[ii] > 1.0){
      // should not happen in the current table
      printf(" *** Warning : extent interpolation becomes extrapolation: ifac_x = %.1f \n", ipol.ifac_x[ii]);
      printf("Better limit a range of values within table limits \n");
    }
  }
}


/**
 * @brief determine the gshift^Gamma for table grid values. Required for further interpolation
 */
void calc_energy_shift_gamma_tabulated(extendedSourceData &dat_ind_a, const int ind_h[2], const int ind_x[2],
                                       double ****energy_shift_bins, int num_bins, int num_shift_bins, double gamma){
  for (int i : {ind_h[0], ind_h[0] + 1, ind_h[1], ind_h[1] + 1}) { // vertical bins, ind_h
    for (int j : {ind_x[0], ind_x[0] + 1, ind_x[1], ind_x[1] + 1}) { // radial extent bins, ind_r
      for (int k = 0; k < num_bins; k++) { // disk radial bins or lensing inclination bins
        double sum = 0.0; // sum is zeroed every time
        for (int l = 0; l < num_shift_bins; l++) { // energy shift bins
          sum += pow(energy_shift_bins[i][j][k][l], gamma);
          assert(sum >= 0.0);
        }
        dat_ind_a.gshift[i][j][k] = sum / num_shift_bins;
      }
    }
  }
}

// a bit of duplicated functions for now. I am not sure what is a better way to avoid duplications
//TODO avoid duplications
double calc_energy_shift_gamma(extendedSourceTable *tab, InterpFactors ipol, double bin_val, double gamma) {
  // do interpolation between bin values due to the fact: if binning at different spins/whatever is different --
  // we have a slight shift of obtained values. In perspective would be good to divide these cases
  int ind_bin_lo, ind_bin_hi;
  double ifac_bin_lo, ifac_bin_hi;
  get_ipol_factor_double(bin_val, tab->dat[ipol.ind_a]->bin, tab->n_bin, &ind_bin_lo, &ifac_bin_lo);
  get_ipol_factor_double(bin_val, tab->dat[ipol.ind_a + 1]->bin, tab->n_bin, &ind_bin_hi, &ifac_bin_hi);

  // some extra operations happen here, if we need only one bin in the end, not all
  //calc_energy_shift_gamma_tabulated(*tab->dat[ipol.ind_a], ipol.ind_h, ipol.ind_x,
                                    //tab->dat[ipol.ind_a]->energy_data_storage,
                                    //tab->n_bin, tab->n_storage, gamma);
  //calc_energy_shift_gamma_tabulated(*tab->dat[ipol.ind_a + 1], ipol.ind_h, ipol.ind_x,
                                    //tab->dat[ipol.ind_a + 1]->energy_data_storage,
                                    //->n_bin, tab->n_storage, gamma);

  for (int a_i : {ipol.ind_a, ipol.ind_a + 1}) { // spin bins
    for (int i : {ipol.ind_h[0], ipol.ind_h[0] + 1, ipol.ind_h[1], ipol.ind_h[1] + 1}) { // vertical bins
      for (int j : {ipol.ind_x[0], ipol.ind_x[0] + 1, ipol.ind_x[1], ipol.ind_x[1] + 1}) { // radial extent bins
        if (a_i == ipol.ind_a) {
          for (int k : {ind_bin_lo, ind_bin_lo + 1}) { // disk radial bins or lensing inclination bins
            double sum = 0.0; // sum is zeroed every time
            for (int l = 0; l < tab->n_storage; l++) { // energy shift bins
              sum += pow(tab->dat[a_i]->energy_data_storage[i][j][k][l], gamma);
              assert(sum >= 0.0);
            }
            tab->dat[a_i]->gshift[i][j][k] = sum / tab->n_storage;
          }
        }
        else {
          for (int k : {ind_bin_hi, ind_bin_hi + 1}) {
            double sum = 0.0;
            for (int l = 0; l < tab->n_storage; l++) {
              sum += pow(tab->dat[a_i]->energy_data_storage[i][j][k][l], gamma);
              assert(sum >= 0.0);
            }
            tab->dat[a_i]->gshift[i][j][k] = sum / tab->n_storage;
          }
        }
      }
    }
  }

  double tab_lo_00 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_bin_lo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_bin_lo + 1]);
  double tab_lo_01 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_bin_lo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_bin_lo + 1]);
  double tab_lo_10 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_bin_lo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_bin_lo + 1]);
  double tab_lo_11 = interp_lin_1d(ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_bin_lo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_bin_lo + 1]);
  double tab_hi_00 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_bin_hi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_bin_hi + 1]);
  double tab_hi_01 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_bin_hi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_bin_hi + 1]);
  double tab_hi_10 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_bin_hi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_bin_hi + 1]);
  double tab_hi_11 = interp_lin_1d(ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_bin_hi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_bin_hi + 1]);
/*
  double g_gamma_low = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_bin],
                                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_bin],
                                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_bin],
                                              tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_bin],
                                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_bin],
                                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_bin],
                                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_bin],
                                              tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_bin]);
  double g_gamma_high = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                               tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_bin + 1],
                                               tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_bin + 1],
                                               tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_bin + 1],
                                               tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_bin + 1],
                                               tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_bin + 1],
                                               tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_bin + 1],
                                               tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_bin + 1],
                                               tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_bin + 1]);

  //double g_gamma = interp_lin_1d(ifac_bin, g_gamma_low, g_gamma_high);
*/
  double g_gamma = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                          tab_lo_00,tab_lo_01,tab_lo_10,tab_lo_11,
                                          tab_hi_00,tab_hi_01,tab_hi_10,tab_hi_11);
  return g_gamma;
}



/**
 * @brief Calculates reflection fraction from interpolated table values, and photon fractions
 * @details
 * @param param, const relParam, a pointer to the structure of parameters used, the parameters include:
 * a - spin, height - height, r_ring - ring radius
 * @param tab, extendedSourceTable, a pointer to the structure containing a table
 * @return lpReflFrac, a pointer to the structure containing interpolated reflection fraction and photon fractions
 */
static lpReflFrac *calc_refl_frac_ext(const extendedSourceTable *tab, InterpFactors ipol, int *status) {

  // Note that this function partially repeats the interpol_estable function. Maybe combine them?
  // the difference is mainly the return.

  lpReflFrac *fractions = new_lpReflFrac(status);
  CHECK_STATUS_RET(*status, fractions);

  fractions->f_ad = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                           tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0]][ipol.ind_x[0]],
                                           tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0]][ipol.ind_x[0] + 1],
                                           tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0] + 1][ipol.ind_x[0]],
                                           tab->dat[ipol.ind_a]->f_ad[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1],
                                           tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1]][ipol.ind_x[1]],
                                           tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1]][ipol.ind_x[1] + 1],
                                           tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1] + 1][ipol.ind_x[1]],
                                           tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1]);

  fractions->f_inf = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                            tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0]][ipol.ind_x[0]],
                                            tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0]][ipol.ind_x[0] + 1],
                                            tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0] + 1][ipol.ind_x[0]],
                                            tab->dat[ipol.ind_a]->f_inf[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1],
                                            tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1]][ipol.ind_x[1]],
                                            tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1]][ipol.ind_x[1] + 1],
                                            tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1] + 1][ipol.ind_x[1]],
                                            tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1]);
  // new table contains this value too
  fractions->f_bh = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                           tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0]][ipol.ind_x[0]],
                                           tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0]][ipol.ind_x[0] + 1],
                                           tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0] + 1][ipol.ind_x[0]],
                                           tab->dat[ipol.ind_a]->f_bh[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1],
                                           tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1]][ipol.ind_x[1]],
                                           tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1]][ipol.ind_x[1] + 1],
                                           tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1] + 1][ipol.ind_x[1]],
                                           tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1]);

  // before there was one more interpolation of refl_frac between values written in function loadSingle source ... above
  // but I don't see any obstacles against just calculating it here. Need to check if it is calculated the same
  fractions->refl_frac = fractions->f_ad / fractions->f_inf;

  // photons are not allowed to cross the disk plane
  if (fractions->f_inf > 0.5) { // should never be the case - TODO in new table
    fractions->f_inf = 0.5;
    printf("Rellp_Extended: inf fraction more than 0.5. \n");
    printf("Replacing with 0.5. This will be fixed in the new extended table. \n");
  }

  ////str->f_bh = 1.0 - (str->f_ad + str->f_inf); // this is not true if we re-define f_inf to be <= 0.5.
  // But it is important only for very large heights.
  // TODO: This definition does not take into account a possible variability of the disk edges. Is there any way to
  // overcome this limitation? Probably only approximate way, because with raytracing it is not easy to vary r_in/r_out


  /// these are equal for now. f_inf_rest is used to calculate primary spectrum
  /// (function flux_boost_source_to_observer in PrimarySource.h)
  /// later maybe would be good to distinguish them correctly - AN
  fractions->f_inf_rest = fractions->f_inf; /// anyway this is true for non-moving sources (beta=0)
  // str->f_inf_rest = 0.5 * (1.0 + cos(del_emit_ad_max));

  return fractions;
}

/*
 * EXTENDED SOURCE GEOMETRY --- MAIN ROUTINE
 * */
/**
 * @brief The main routine
 * @details This routine calculates emissivity for the extended source using calc_extended_emis function declared further
 * @param emis, emisProfile, pointer to the structure containing emissivity profile
 * @param param, relParam, pointer to the structure containing the parameters used
 */
void get_emis_extendedSource(emisProfile *emis, const relParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  extendedSourceTable* emissivity_table = get_es_table(status);
  // new part - loading primary normalization, i.e., lensing_factor * g_{so}^Gamma, to primary normalization
  extendedSourceTable* prim_table = get_prim_table(status);
  calc_extended_emis(emis, param, emissivity_table, prim_table, status); // PRIM , prim_table

}

/**
 * @brief Secondary main routine which gets g^Gamma from source to observer from the table (at relxill inclination)
 * @details for now it is almost duplicated with the function below, which does the same for the disk point
 * @return The g^Gamma value with Gamma not necessarily taken from relParam (e.g. can be 1 for just energyshift)
 */
double get_energyboost_ext_source_obs(const relParam *param, double gamma) {
  int status = EXIT_SUCCESS;
  extendedSourceTable* prim_table = get_prim_table(&status); // important to load the correct table
  InterpFactors prim_ipol; // just get the factors to the closest table grid points in a, h, x dimensions
  get_ipol_factors_a_h_x(param, prim_table, prim_ipol);
  // input incl can be both in radian and in grads. The table has cos incl, so we (non-optimally) convert it here
  // (I use the fact that relxill limits are 3 and 87 deg, so if input incl is < 1.57 (~pi/2) - it is radian)
  double cos_incl = (param->incl < 0.5 * M_PI) ? cos(param->incl) : cos(param->incl * M_PI / 180.0);
  // in the first case it is in radian, in the second it is in grad
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
  int status = EXIT_SUCCESS; // dumb solution to re-define status here, but we don't have it in the Relphysics functions
  // and we have it in getting table routines...
  extendedSourceTable* table = get_es_table(&status); // important to load the correct table
  InterpFactors ipol; // just get the factors to the closest table grid points in a, h, x dimensions
  get_ipol_factors_a_h_x(param, table, ipol);
  double gshift = calc_energy_shift_gamma(table, ipol, radius, gamma);
  return gshift;
}


/**
 * @brief reads extended source table
 * @details
 * @param filename; const char, pointer to the name of the file (now declared in Rellp_Extended.h)
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
      tab = new_extendedSourceTable(PRIM_TABLE_NHEIGHT,PRIM_TABLE_NDISTANCE, PRIM_TABLE_NA,
                                    PRIM_TABLE_NINCL, PRIM_TABLE_NENERGYBINS, status);
    } else {
      tab = new_extendedSourceTable(EXT_TABLE_NHEIGHT, EXT_TABLE_NDISTANCE, EXT_TABLE_NA,
                                    EXT_TABLE_NRAD, EXT_TABLE_NENERGYBINS, status);
    }

    CHECK_STATUS_BREAK(*status);

    // should be set by previous routine
    assert(tab != nullptr);
    assert(tab->dat != nullptr);

    // get the full filename
    if (sprintf(fullfilename, "%s/%s", get_relxill_table_path(), filename) == -1) {
      RELXILL_ERROR("failed to construct full path the rel table\n", status);
      break;
    }

    // open the file
    if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
     // fits_report_error(stdout, *status);
      CHECK_RELXILL_ERROR("opening of the extended source table failed", status);
      printf("    full path given: %s \n", fullfilename);
      break;
    }

    // first read the spin axes of the table AND move to the correct extension
    // printf("************************\n");
    // printf("filename : %s\n", fullfilename);


    get_extendedsourcetable_spin(tab->n_a, &(tab->a) , fptr, status);
    CHECK_RELXILL_ERROR("reading of spin axis failed\n", status);

    // now load the full table (need to go through all extensions)
    int ii;
    int rownum;
    for (ii = 0; ii < tab->n_a; ii++) {

      rownum = ii + 2; // number of the row we read in the fits table
      fits_movabs_hdu(fptr, rownum, nullptr, status);

      assert(tab->dat[ii] == nullptr);

      tab->dat[ii] = load_single_extendedSourceData(
          fptr, tab->n_h, tab->n_r, tab->n_bin, tab->n_storage, status,
          tab->a[ii], is_lensing_table);
      if (*status != EXIT_SUCCESS) {
        RELXILL_ERROR("failed to load data from the extended source table into memory\n", status);

        break;
      }
    }

  } while (false);

  if (*status == EXIT_SUCCESS) {
    // assign the value
    (*inp_tab) = tab;
  } else {
    free_extendedSourceTable(tab, is_lensing_table);
  }

  if(fptr != nullptr) {fits_close_file(fptr, status);}

}

double calc_primary_lensing(double incl, extendedSourceTable *prim_tab, InterpFactors ipol) {
  int ind_incl; // get indexes of interpolation
  double ifac_incl;
  get_ipol_factor_double(incl, prim_tab->dat[ipol.ind_a]->bin, prim_tab->n_bin,
                         &ind_incl, &ifac_incl);
  if (ifac_incl < 0.0 or ifac_incl > 1.0){
    // should not happen in the current table
    printf(" *** Warning : inclination interpolation becomes extrapolation: ifac_i = %.1f \n", ifac_incl);
    printf("Better limit a range of values within table limits, 0.025 < cos(incl) < 0.975 \n");
  }

  double lens_low = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                         prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][ind_incl],
                         prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_incl],
                         prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_incl],
                         prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_incl],
                         prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1]][ind_incl],
                         prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_incl],
                         prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_incl],
                         prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_incl]);
  double lens_high = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                           prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][ind_incl + 1],
                                           prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_incl + 1],
                                           prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_incl + 1],
                                           prim_tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_incl + 1],
                                           prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1]][ind_incl + 1],
                                           prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_incl + 1],
                                           prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_incl + 1],
                                           prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_incl + 1]);
  double lensing = interp_lin_1d(ifac_incl, lens_low, lens_high);
  return lensing;
}

/////////////////////////////////////////////////////////
/**
 * @brief Calculates emissivity
 * @details
 * @param emisProf; emisProfile, pointer to the structure with emissivity profile
 * @param param; relParam, pointer to the parameters structure
 * @param tab; extendedSourceTable, pointer to the table structure
 * @param prim_tab; extendedSourceTable, pointer to the primary table
 */
void calc_extended_emis(emisProfile *emisProf, const relParam *param,
                        extendedSourceTable *tab, extendedSourceTable *prim_tab, int *status) {

  CHECK_STATUS_VOID(*status);

  InterpFactors ipol; // new structure with interpolation factors instead of a ton of indexes and factors
  get_ipol_factors_a_h_x(param, tab, ipol);


  // rename function
  emisProfile *emis_profile_table = interpol_estable(tab, ipol, param->gamma, status);
   rebin_emisprofile_on_radial_grid(emisProf, emis_profile_table, status);

  free(emis_profile_table->re); // is not freed by free_emisProfile (Rellp.cpp)
  free_emisProfile(emis_profile_table);

  InterpFactors prim_ipol;
  get_ipol_factors_a_h_x(param, prim_tab, prim_ipol);
  emisProf->photon_fate_fractions = calc_refl_frac_ext(prim_tab, prim_ipol, status);


  // input incl can be both in radian and in grads. The table has cos incl, so we unoptimally convert it here
  // (I use the fact that relxill limits are 3 and 87 deg, so if input incl is < 1.57 (~pi/2) - it is radian)
  double cos_incl = (param->incl < 0.5 * M_PI) ? cos(param->incl) : cos(param->incl * M_PI / 180.0);
  // in the first case it is in radian, in the second it is in grad

    // experimentally put this calculation to refl_frac structure
  emisProf->photon_fate_fractions->lensing = calc_primary_lensing(cos_incl, prim_tab, prim_ipol);

  // Lensing is currently calculated with raytracing for disk size r_out = 10^3 rg.
  // Therefore at large heights we loose some fraction of photons (they go "under the disk")
  // That is why I add here this empirical correction, which plays some role at large h
  emisProf->photon_fate_fractions->lensing /= (1.0 -
      sqrt(param->height * param->height + param->r_ring * param->r_ring + param->a * param->a) / 1000.0);

  // The last step is to add boost g^Gamma to observer, which could also be done in PrimarySource.h
  emisProf->photon_fate_fractions->lensing_and_boost_factor = (emisProf->photon_fate_fractions->lensing *
      calc_energy_shift_gamma(prim_tab, prim_ipol, cos_incl, param->gamma));


  // take the beaming of the jet into account (see Dauser et al., 2013)
  // Not done for extended source, todo
//  for (int ii=0; ii<tab->n_rad; ii++){
//    if (param->beta > 1e-6) {
//      emisProf->emis[ii] *= pow(doppler_factor(emisProf->del_emit[ii], param->beta), 2);
//    }
//  }
}

/**
 * @brief Interpolation of the table between spins
 * @details
 * @param a; double, spin
 * @param height; double, height
 * @param r_ring; double, ring radius
 * @param gamma; double, photon index
 * @param tab; extendedSourceTable, pointer to the table as declared above
 * @return emis_profile_table; emisProfile, pointer to the structure with emissivity profile
 */
emisProfile *interpol_estable(extendedSourceTable *tab, InterpFactors ipol, double gamma, int *status) {
  CHECK_STATUS_RET(*status, nullptr);

  // Need to interpolate radial grid of table because it is different for different spins (risco(a))
  auto* radial_grid = (double*) malloc(sizeof(double)*tab->n_bin);
  CHECK_MALLOC_RET_STATUS(radial_grid,status,nullptr)
  for (int ii = 0; ii < tab->n_bin; ii++) {
    radial_grid[ii] = interp_lin_1d(ipol.ifac_a, tab->dat[ipol.ind_a]->bin[ii], tab->dat[ipol.ind_a+1]->bin[ii]);
  }

  emisProfile* emis_profile_table = new_emisProfile(radial_grid, tab->n_bin, status);

  // It is duplicated and could be done inside of existing element-vise function for g_gamma
  calc_energy_shift_gamma_tabulated(*tab->dat[ipol.ind_a], ipol.ind_h, ipol.ind_x,
                                    tab->dat[ipol.ind_a]->energy_data_storage,
                          tab->n_bin, tab->n_storage, gamma);
  calc_energy_shift_gamma_tabulated(*tab->dat[ipol.ind_a + 1], ipol.ind_h, ipol.ind_x,
                                    tab->dat[ipol.ind_a + 1]->energy_data_storage,
                          tab->n_bin, tab->n_storage, gamma);

  // pre-allocation of values shifted in radial grid before interpolation
  double tab_lo_00, tab_lo_01, tab_lo_10, tab_lo_11, tab_hi_00, tab_hi_01, tab_hi_10, tab_hi_11;
  double fac_rlo, fac_rhi;
  int ind_rlo, ind_rhi;
  for (int ii = 0; ii < tab->n_bin; ii++) {
    // as radial grids shift between spins, we try to compensate for this shift
    // I know that sometimes we can extrapolate if r < r_low_spin. But normally only a bit
    // If extrapolation comes too far -- then this grid is far from the desired spin
    get_ipol_factor_double(radial_grid[ii], tab->dat[ipol.ind_a]->bin, tab->n_bin, &ind_rlo, &fac_rlo);
    get_ipol_factor_double(radial_grid[ii], tab->dat[ipol.ind_a + 1]->bin, tab->n_bin, &ind_rhi, &fac_rhi);
    // Interpolation of the flux
    // looks like meh, but gives some generalization potential (maybe refactor to smth better reading further?)
    //if (fac_rlo < 0.0 or fac_rlo > 1.0) {
      //printf("Warning, radius interpolation out of bounds, at low spin interp factor is: %.6f at radius: %.6f\n", fac_rlo, radial_grid[ii]);
    //}
    //if (fac_rhi < 0.0 or fac_rhi > 1.0) { // shouldn't be for a > 0
      //printf("Warning, radius interpolation out of bounds, at high spin interp factor is: %.6f at radius: %.6f\n", fac_rhi, radial_grid[ii]);
    //}

    tab_lo_00 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo], tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_01 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo], tab->dat[ipol.ind_a]->val[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo + 1]);
    tab_lo_10 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo], tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_11 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo], tab->dat[ipol.ind_a]->val[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo + 1]);
    tab_hi_00 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi], tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_01 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi], tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi + 1]);
    tab_hi_10 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi], tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_11 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi], tab->dat[ipol.ind_a + 1]->val[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi + 1]);
    emis_profile_table->emis[ii] = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                                          tab_lo_00,tab_lo_01,
                                                          tab_lo_10,tab_lo_11,
                                                          tab_hi_00,tab_hi_01,
                                                          tab_hi_10,tab_hi_11);

    // is it safe to reuse these temporary variables twice in cycle?
    tab_lo_00 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_01 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0]][ipol.ind_x[0] + 1][ind_rlo + 1]);
    tab_lo_10 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0]][ind_rlo + 1]);
    tab_lo_11 = interp_lin_1d(fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo], tab->dat[ipol.ind_a]->gshift[ipol.ind_h[0] + 1][ipol.ind_x[0] + 1][ind_rlo + 1]);
    tab_hi_00 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_01 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1]][ipol.ind_x[1] + 1][ind_rhi + 1]);
    tab_hi_10 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1]][ind_rhi + 1]);
    tab_hi_11 = interp_lin_1d(fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi], tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_h[1] + 1][ipol.ind_x[1] + 1][ind_rhi + 1]);

    double g_gamma = get_interpolated_value(ipol.ifac_a, ipol.ifac_h, ipol.ifac_x,
                                            tab_lo_00,tab_lo_01,
                                            tab_lo_10,tab_lo_11,
                                            tab_hi_00,tab_hi_01,
                                            tab_hi_10,tab_hi_11);


    // this line accounts for g^Gamma:
    emis_profile_table->emis[ii] *= g_gamma; // in fact, this is "boost" of the flux, which for lp model is taken
    // into account only in calc_fluxboost_source_disk
    assert(emis_profile_table->emis[ii] >= 0.0); // if interpolation behaves wrong
  }
  return emis_profile_table;
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
  tab->n_r = n_r;

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
extendedSourceData *new_extendedSourceData(int *status, bool lensing) {
  int n_h = EXT_TABLE_NHEIGHT;
  int n_x = EXT_TABLE_NDISTANCE;
  int n_bin = EXT_TABLE_NRAD;
  int n_shift = EXT_TABLE_NENERGYBINS;
  if (lensing) {
    n_h = PRIM_TABLE_NHEIGHT;
    n_x = PRIM_TABLE_NDISTANCE;
    n_bin = PRIM_TABLE_NINCL;
    n_shift = PRIM_TABLE_NENERGYBINS;
  }
  auto *dat = (extendedSourceData *) malloc(sizeof(extendedSourceData));
  CHECK_MALLOC_RET_STATUS(dat, status, nullptr)

  dat->height = (double *) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->height, status, nullptr)

  dat->r_ring = (double **) malloc(sizeof(double) * n_h);
  CHECK_MALLOC_RET_STATUS(dat->r_ring, status, nullptr)

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

  for (int i=0; i<EXT_TABLE_NHEIGHT; i++){

    dat->refl_frac[i] = (double *) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, nullptr)

    dat->r_ring[i] = (double *) malloc(sizeof(double) * n_x);
    CHECK_MALLOC_RET_STATUS(dat->r_ring, status, nullptr)

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
 *
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
 *
 * @brief Frees the memory stored before for the data
 * @details
 * @param dat; extendedSourceData, data to be freed
 */
void free_extendedSourceData(extendedSourceData *dat, bool is_lensing_table) {
  int n_h = EXT_TABLE_NHEIGHT;
  int n_x = EXT_TABLE_NDISTANCE;
  int n_r = EXT_TABLE_NRAD;
  if (is_lensing_table) {
    n_h = PRIM_TABLE_NHEIGHT; // re-define in case of lensing table, dimensions can be different
    n_x = PRIM_TABLE_NDISTANCE;
    n_r = PRIM_TABLE_NINCL;
  }
  if (dat != nullptr) {
    for (int ii=0; ii < n_h; ii++) {
      if (dat->r_ring != nullptr) free(dat->r_ring[ii]);
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
    free(dat->r_ring);
    free(dat->val);
    free(dat->gshift);
    free(dat->energy_data_storage);
    free(dat->refl_frac);
    free(dat->f_inf);
    free(dat->f_ad);
    free(dat->f_bh);

}
