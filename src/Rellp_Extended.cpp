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

extendedSourceTable* cached_extendedSource_table = nullptr; // emissivity table

extendedSourceTable* cached_extendedSource_table_prim = nullptr; // lensing table

/**
 * @brief Gets extended source table with a given name EXT_TABLE_FILENAME
 * @details
 * @param status pointer, status of the code for debugging
 * @return cached_extendedSource_table
 */
extendedSourceTable* get_es_table(int* status) {
  CHECK_STATUS_RET(*status, nullptr);

  if (cached_extendedSource_table == nullptr) {
    read_extended_source_table(EXT_TABLE_FILENAME, &cached_extendedSource_table, status);
    CHECK_STATUS_RET(*status, nullptr);
  }
  return cached_extendedSource_table;
}

extendedSourceTable* get_prim_table(int* status) {
  CHECK_STATUS_RET(*status, nullptr);

  bool is_lensing_table = true; // quick solution to have two tables instead of one
  if (cached_extendedSource_table_prim == nullptr) {
    read_extended_source_table(PRIM_TABLE_FILENAME, &cached_extendedSource_table_prim, status,
                               is_lensing_table);
    CHECK_STATUS_RET(*status, nullptr);
  }
  return cached_extendedSource_table_prim;
}

/**
 * @brief Gets rows of the table corresponding to different spins
 * @details Uses cfitsio routines to read preloaded fits file
 * @param n_rows integer, number of rows in the table
 * @param val two-fold pointer to the value which stores the table rows
 * @param fptr pointer to the file
 * @param status pointer status of the code for debugging
 */
static void get_spin_rows_extended_table(const int n_rows, float** val, fitsfile* fptr, int* status) {

  *val = (float*)  malloc(n_rows * sizeof(float));

  for (int i_r = 0; i_r < n_rows; i_r++){

    int index = i_r + 2; // number of the row we read in the fits table
    fits_movabs_hdu(fptr, index, nullptr, status);

    fits_read_key(fptr, TFLOAT, "SPIN", &((*val)[i_r]), nullptr, status);
    if (*status != EXIT_SUCCESS) {
      CHECK_RELXILL_ERROR("reading spin dimension in table failed,", status);
    }
  }

}


/**
 * @brief Gets data for a single source from the table and then makes a relativistic correction of the flux
 * @details
 * @param fptr pointer, to the fits file location
 * @param n_r integer, number of spherical radii simulated
 * @param n_th integer, number of polar angles simulated
 * @param n_bin integer, number of radial bins
 * @param n_storage integer, number of energyshift storage bins. Used only to read the fifth column with shift
 * @param spin double, spin value of the source
 * @return extendedSourceData, data from the table
 */
static extendedSourceData* load_single_extended_source_data(
    fitsfile* fptr, const int n_r, const int n_th, const int n_bin, const int n_storage,
    int* status, const double spin, const bool is_lensing_table) {

  extendedSourceData *dat = new_extended_source_data(status, is_lensing_table);
  CHECK_MALLOC_RET_STATUS(dat, status, nullptr)

  int anynullptr = 0;
  double doublenullptr = 0.0;

  // read the columns
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, n_bin, &doublenullptr,
                dat->bin,&anynullptr, status);

  for (int i_r = 0; i_r < n_r; i_r++) {
    fits_read_col(fptr, TDOUBLE, 2, 1 + i_r * n_th, 1, 1, &doublenullptr,
                  &(dat->r_sph[i_r]), &anynullptr, status);
    fits_read_col(fptr, TDOUBLE, 3, 1 + i_r * n_th, 1, n_th,
                  &doublenullptr, dat->theta[i_r], &anynullptr, status);
  }

  for (int i_r = 0; i_r < n_r; i_r++) {
    for (int j_th = 0; j_th < n_th; j_th++) {

      int i_row = (i_r * n_th) + j_th + 1;

      fits_read_col(fptr, TDOUBLE, 4, i_row, 1, n_bin, &doublenullptr,
                    dat->val[i_r][j_th], &anynullptr, status);

      for (int k_b = 0; k_b < n_bin; k_b++) {
        // should never happen with the proper table:
        if (dat->val[i_r][j_th][k_b] < 0.0) {
          CHECK_RELXILL_ERROR("Negative value loaded from table, check your table", status);
        }
        fits_read_col(fptr, TDOUBLE, 5, i_row, k_b * n_storage + 1, n_storage,
                      &doublenullptr, dat->energy_data_storage[i_r][j_th][k_b], &anynullptr, status);
      }

      // n_rad to 1 because it is not an array just a single value
      fits_read_col(fptr, TDOUBLE, 6, i_row, 1, 1, &doublenullptr,
                    &dat->f_ad[i_r][j_th], &anynullptr, status);

      fits_read_col(fptr, TDOUBLE, 7, i_row, 1, 1, &doublenullptr,
                    &dat->f_inf[i_r][j_th], &anynullptr, status);
      if (is_lensing_table) { // extra column which is now only in the primary table
        fits_read_col(fptr, TDOUBLE, 8, i_row, 1, 1, &doublenullptr,
                      &dat->f_bh[i_r][j_th], &anynullptr, status);
      }
      dat->refl_frac[i_r][j_th] = dat->f_ad[i_r][j_th] / dat->f_inf[i_r][j_th];

    }
  }
  if (!is_lensing_table) { // only for flux, not for lensing_factor
    apply_flux_correction(dat, dat->bin, spin, n_r, n_th, n_bin);
  }
  return dat;
}


/**
 * @brief Reads extended source table
 * @details
 * @param filename; const char, pointer to the name of the file (declared in Rellp_Extended.h)
 * @param inp_tab; extendedSourceTable, two-fold pointer to the input table data
 */
void read_extended_source_table(const char* filename, extendedSourceTable** inp_tab, int* status,
                                const bool is_lensing_table) {

  extendedSourceTable *tab = (*inp_tab);
  fitsfile *fptr = nullptr;

  char fullfilename[999];

  do {
    if (tab != nullptr) {
      RELXILL_ERROR("extended source table already loaded\n", status);
      break;
    }
    if (is_lensing_table) {
      tab = new_extended_source_table(PRIM_TABLE_NUM_RADII, PRIM_TABLE_NUM_ANGLES, PRIM_TABLE_NUM_SPINS,
                                      PRIM_TABLE_NUM_INCLINATIONS, PRIM_TABLE_NUM_ENERGY_BINS, status);
    } else {
      tab = new_extended_source_table(EXT_TABLE_NUM_RADII, EXT_TABLE_NUM_ANGLES, EXT_TABLE_NUM_SPINS,
                                      EXT_TABLE_NUM_RAD_BINS, EXT_TABLE_NUM_ENERGY_BINS, status);
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
    int i_row;
    for (int i_a = 0; i_a < tab->n_a; i_a++) {
      i_row = i_a + 2; // number of the row we read in the fits table
      fits_movabs_hdu(fptr, i_row, nullptr, status);

      if (tab->dat[i_a] != nullptr) {
        CHECK_RELXILL_ERROR("reading of spin axis failed\n", status);
      }

      tab->dat[i_a] = load_single_extended_source_data(
          fptr, tab->n_r, tab->n_th, tab->n_bin, tab->n_storage, status, tab->a[i_a], is_lensing_table);
      if (*status != EXIT_SUCCESS) {
        RELXILL_ERROR("failed to load data from the extended source table into memory\n", status);

        break;
      }
    }

  } while (false);

  if (*status == EXIT_SUCCESS) {
    (*inp_tab) = tab;
  } else {
    free_extended_source_table(tab, is_lensing_table);
  }

  if(fptr != nullptr) {fits_close_file(fptr, status);}
}


/**
 * @brief Interpolate the table between two spins
 * @details
 * @param tab; extendedSourceTable, pointer to the table as declared above
 * @param ipol; structure with interpolation factors between a, r, theta
 * @param gamma; double, photon index
 * @return emis_profile_table; emisProfile, pointer to the structure with emissivity profile
 */
emisProfile* interpol_extended_table(extendedSourceTable* tab, InterpFactors ipol,
                                     const double gamma, double* radial_grid, int* status) {
  CHECK_STATUS_RET(*status, nullptr);

  // The radial grid of the table is different for different spins (risco(a)), so we must interpolate it
  assert(radial_grid == nullptr);
  radial_grid = (double*) malloc(sizeof(double)*tab->n_bin);
  CHECK_MALLOC_RET_STATUS(radial_grid,status,nullptr)
  for (int ii = 0; ii < tab->n_bin; ii++) {
    radial_grid[ii] = interp_lin_1d(ipol.ifac_a,
                                    tab->dat[ipol.ind_a]->bin[ii],tab->dat[ipol.ind_a+1]->bin[ii]);
  }

  emisProfile* emis_profile_table = new_emisProfile(radial_grid, tab->n_bin, status);

  calc_energy_shift_gamma_tabulated(
      *tab->dat[ipol.ind_a], ipol.ind_r[0], ipol.ind_th[0],
      tab->dat[ipol.ind_a]->energy_data_storage,tab->n_bin, tab->n_storage,
      gamma);
  calc_energy_shift_gamma_tabulated(
      *tab->dat[ipol.ind_a + 1], ipol.ind_r[1], ipol.ind_th[1],
      tab->dat[ipol.ind_a + 1]->energy_data_storage,tab->n_bin, tab->n_storage,
      gamma);

  // pre-allocation of values shifted in radial grid before interpolation
  double tab_lo_00, tab_lo_01, tab_lo_10, tab_lo_11, tab_hi_00, tab_hi_01, tab_hi_10, tab_hi_11;
  double fac_rlo, fac_rhi;
  int ind_rlo, ind_rhi;
  for (int ii = 0; ii < tab->n_bin; ii++) {
    // As radial grids shift between spins, we try to compensate for this shift in the values
    get_ipol_factor_double(radial_grid[ii], tab->dat[ipol.ind_a]->bin,
                           tab->n_bin, &ind_rlo, &fac_rlo);
    get_ipol_factor_double(radial_grid[ii], tab->dat[ipol.ind_a + 1]->bin,
                           tab->n_bin, &ind_rhi, &fac_rhi);

    // Interpolation of the flux
    tab_lo_00 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0]][ind_rlo],
        tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0]][ind_rlo + 1]);
    tab_lo_01 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_rlo],
        tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_rlo + 1]);
    tab_lo_10 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_rlo],
        tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_rlo + 1]);
    tab_lo_11 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_rlo],
        tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_rlo + 1]);

    tab_hi_00 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1]][ind_rhi],
        tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1]][ind_rhi + 1]);
    tab_hi_01 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_rhi],
        tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_rhi + 1]);
    tab_hi_10 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_rhi],
        tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_rhi + 1]);
    tab_hi_11 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_rhi],
        tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_rhi + 1]);

    emis_profile_table->emis[ii] = get_interpolated_value_a_r_th(
        ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
        tab_lo_00, tab_lo_01,tab_lo_10, tab_lo_11,
        tab_hi_00, tab_hi_01,tab_hi_10, tab_hi_11);

    // Re-use temporary variables for both emis and g_gamma
    tab_lo_00 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0]][ind_rlo],
        tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0]][ind_rlo + 1]);
    tab_lo_01 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_rlo],
        tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_rlo + 1]);
    tab_lo_10 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_rlo],
        tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_rlo + 1]);
    tab_lo_11 = interp_lin_1d(
        fac_rlo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_rlo],
        tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_rlo + 1]);

    tab_hi_00 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1]][ind_rhi],
        tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1]][ind_rhi + 1]);
    tab_hi_01 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_rhi],
        tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_rhi + 1]);
    tab_hi_10 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_rhi],
        tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_rhi + 1]);
    tab_hi_11 = interp_lin_1d(
        fac_rhi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_rhi],
        tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_rhi + 1]);

    double g_gamma = get_interpolated_value_a_r_th(
        ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
        tab_lo_00, tab_lo_01,tab_lo_10, tab_lo_11,
        tab_hi_00, tab_hi_01,tab_hi_10, tab_hi_11);

    emis_profile_table->emis[ii] *= g_gamma;
    if (emis_profile_table->emis[ii] < 0.0) {
      CHECK_RELXILL_ERROR("Negative emissivity bin caught after table interpolation!", status);
    }
  }
  return emis_profile_table;
}


/**
 * @brief Interpolates a given tabulated value (between 4 base points) between spherical radii and polar angles,
 * and then between spins. In fact, this is 3d interpolation in (a,r,th)
 * @return interpolated value, double
 */
double get_interpolated_value_a_r_th(
    const double ifac_a, const double ifac_r[2], const double ifac_th[2],
    const double vals_lo00, const double vals_lo01, const double vals_lo10, const double vals_lo11,
    const double vals_hi00, const double vals_hi01, const double vals_hi10, const double vals_hi11) {
  // does h-r interpolation for spin[ind_a]
  double val_at_spin_lo = interp_lin_2d(ifac_th[0], ifac_r[0],
                                        vals_lo00, vals_lo01, vals_lo10, vals_lo11);
  // does h-r interpolation for spin[ind_a+1]
  double val_at_spin_hi = interp_lin_2d(ifac_th[1], ifac_r[1],
                                        vals_hi00, vals_hi01, vals_hi10, vals_hi11);
  // grid points must be equal, no adjustment done here
  return interp_lin_1d(ifac_a, val_at_spin_lo, val_at_spin_hi);
}


void get_ipol_factors_a_r_th(const double a, const double r, const double th,
                             const extendedSourceTable* tab, InterpFactors& ipol) {
  // get interpolation factors between spins
  get_ipol_factor((float) a, tab->a, tab->n_a, &ipol.ind_a, &ipol.ifac_a);
  if (ipol.ifac_a < 0.0 or ipol.ifac_a > 1.0) {
  }
  for (int ii = 0; ii < 2; ii++) {
    // get interpolation factors between spherical radii
    get_ipol_factor_double(r, tab->dat[ipol.ind_a + ii]->r_sph, tab->n_r,
                           &(ipol.ind_r[ii]), &(ipol.ifac_r[ii]));
    if (ipol.ifac_r[ii] < 0.0 or ipol.ifac_r[ii] > 1.0) {
      printf(" *** Warning : spherical radius interpolation becomes extrapolation: ifac_r = %.1f \n",
             ipol.ifac_r[ii]);
      printf(" *** The values: a = %.3f, r_src = %.3f, theta_src = %.3f \n", a, r, th);
      printf("Better limit a range of values within table limits \n");
    }
    // interpolation factors between polar angles
    int index = (ipol.ifac_r[ii] > 0.999999) ? ipol.ind_r[ii] + 1 : ipol.ind_r[ii];
    // this additional index shift doesn't matter for the current table
    // but if tabulated angles will differ for different radii - it would matter
    get_ipol_factor_double(th, tab->dat[ipol.ind_a + ii]->theta[index], tab->n_th,
                           &(ipol.ind_th[ii]), &(ipol.ifac_th[ii]));
    if (ipol.ifac_th[ii] < 0.0 or ipol.ifac_th[ii] > 1.0) {
      printf(" *** Warning : angle interpolation becomes extrapolation: ifac_th = %.1f \n", ipol.ifac_th[ii]);
      printf("Better limit a range of values within table limits \n");
    }
  }
}


/**
 * @brief determine the gshift^Gamma for table grid values. Required for further interpolation
 */
void calc_energy_shift_gamma_tabulated(extendedSourceData& dat_ind_a, const int ind_r, const int ind_th,
                                       double**** energy_shift_bins, const int num_bins, const int num_shift_bins,
                                       const double gamma) {
  for (int i_r : {ind_r, ind_r + 1}) { // sph radial bins
    for (int j_th : {ind_th, ind_th + 1}) { // polar angle bins
      for (int k_b = 0; k_b < num_bins; k_b++) { // disk radial bins or lensing inclination bins
        double _sum = 0.0; // sum is zeroed every time
        for (int l_e = 0; l_e < num_shift_bins; l_e++) { // energy shift bins
          _sum += pow(energy_shift_bins[i_r][j_th][k_b][l_e], gamma);
          assert(_sum >= 0.0 || !"summation of g is negative, wrong value in table\n");
        }
        dat_ind_a.gshift[i_r][j_th][k_b] = _sum / num_shift_bins;
      }
    }
  }
}


/**
 * @brief determine the gshift^Gamma for an arbitrary point (a, r, theta)
 * @details The difference to calc_energy_shift_gamma_tabulated function is that here we are not on the grid points
 * of the table, so we need to do extra actions. This function is more general than calc_energy_shift_gamma_tabulated
 * But obviously a bit slower too due to extra actions
 */
double calc_energy_shift_gamma(
    extendedSourceTable* tab, const InterpFactors ipol, const double bin_val, const double gamma) {
  int ind_bin_lo, ind_bin_hi;
  double ifac_bin_lo, ifac_bin_hi;
  get_ipol_factor_double(bin_val, tab->dat[ipol.ind_a]->bin,
                         tab->n_bin, &ind_bin_lo, &ifac_bin_lo);
  get_ipol_factor_double(bin_val, tab->dat[ipol.ind_a + 1]->bin,
                         tab->n_bin, &ind_bin_hi, &ifac_bin_hi);

  for (int i_r : {ipol.ind_r[0], ipol.ind_r[0] + 1}) {
    for (int j_th : {ipol.ind_th[0], ipol.ind_th[0] + 1}) {
      for (int k_b : {ind_bin_lo, ind_bin_lo + 1}) {
        double _sum = 0.0;
        for (int l_e = 0; l_e < tab->n_storage; l_e++) {
          _sum += pow(tab->dat[ipol.ind_a]->energy_data_storage[i_r][j_th][k_b][l_e], gamma);
          assert(_sum >= 0.0 || !"summation of g is negative, wrong value in table\n");
        }
        tab->dat[ipol.ind_a]->gshift[i_r][j_th][k_b] = _sum / tab->n_storage;
      }
    }
  }
  for (int i_r : {ipol.ind_r[1], ipol.ind_r[1] + 1}) {
    for (int j_th : {ipol.ind_th[1], ipol.ind_th[1] + 1}) {
      for (int k_b : {ind_bin_hi, ind_bin_hi + 1}) {
        double _sum = 0.0;
        for (int l_e = 0; l_e < tab->n_storage; l_e++) {
          _sum += pow(tab->dat[ipol.ind_a + 1]->energy_data_storage[i_r][j_th][k_b][l_e], gamma);
          assert(_sum >= 0.0 || !"summation of g is negative, wrong value in table\n");
        }
        tab->dat[ipol.ind_a + 1]->gshift[i_r][j_th][k_b] = _sum / tab->n_storage;
      }
    }
  }

  double tab_lo_00 = interp_lin_1d(
      ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0]][ind_bin_lo],
      tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0]][ind_bin_lo + 1]);
  double tab_lo_01 = interp_lin_1d(
      ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_bin_lo],
      tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_bin_lo + 1]);
  double tab_lo_10 = interp_lin_1d(
      ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_bin_lo],
      tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_bin_lo + 1]);
  double tab_lo_11 = interp_lin_1d(
      ifac_bin_lo, tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_bin_lo],
      tab->dat[ipol.ind_a]->gshift[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_bin_lo + 1]);

  double tab_hi_00 = interp_lin_1d(
      ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1]][ind_bin_hi],
      tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1]][ind_bin_hi + 1]);
  double tab_hi_01 = interp_lin_1d(
      ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_bin_hi],
      tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_bin_hi + 1]);
  double tab_hi_10 = interp_lin_1d(
      ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_bin_hi],
      tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_bin_hi + 1]);
  double tab_hi_11 = interp_lin_1d(
      ifac_bin_hi, tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_bin_hi],
      tab->dat[ipol.ind_a + 1]->gshift[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_bin_hi + 1]);

  double g_gamma = get_interpolated_value_a_r_th(
      ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
      tab_lo_00, tab_lo_01, tab_lo_10, tab_lo_11,
      tab_hi_00, tab_hi_01, tab_hi_10, tab_hi_11);
  assert (g_gamma >= 0.0 || !"g^Gamma negative, wrong interpolation!\n");
  return g_gamma;
}


/**
 * @brief Calculates reflection fraction from interpolated table values, and photon fractions
 * @details Can it be combined/restructured with interpol_extended_table ?
 * @return lpReflFrac, pointer to the structure containing interpolated reflection fraction and photon fractions
 */
static lpReflFrac* calc_refl_frac_ext(const extendedSourceTable* tab, const InterpFactors ipol, int* status) {

  lpReflFrac *fractions = new_lpReflFrac(status);
  CHECK_STATUS_RET(*status, fractions);

  fractions->f_ad = get_interpolated_value_a_r_th(
      ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
      tab->dat[ipol.ind_a]->f_ad[ipol.ind_r[0]][ipol.ind_th[0]],
      tab->dat[ipol.ind_a]->f_ad[ipol.ind_r[0]][ipol.ind_th[0] + 1],
      tab->dat[ipol.ind_a]->f_ad[ipol.ind_r[0] + 1][ipol.ind_th[0]],
      tab->dat[ipol.ind_a]->f_ad[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1],
      tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_r[1]][ipol.ind_th[1]],
      tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_r[1]][ipol.ind_th[1] + 1],
      tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_r[1] + 1][ipol.ind_th[1]],
      tab->dat[ipol.ind_a + 1]->f_ad[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1]);

  fractions->f_inf = get_interpolated_value_a_r_th(
      ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
      tab->dat[ipol.ind_a]->f_inf[ipol.ind_r[0]][ipol.ind_th[0]],
      tab->dat[ipol.ind_a]->f_inf[ipol.ind_r[0]][ipol.ind_th[0] + 1],
      tab->dat[ipol.ind_a]->f_inf[ipol.ind_r[0] + 1][ipol.ind_th[0]],
      tab->dat[ipol.ind_a]->f_inf[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1],
      tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_r[1]][ipol.ind_th[1]],
      tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_r[1]][ipol.ind_th[1] + 1],
      tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_r[1] + 1][ipol.ind_th[1]],
      tab->dat[ipol.ind_a + 1]->f_inf[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1]);

  fractions->f_bh = get_interpolated_value_a_r_th(
      ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
      tab->dat[ipol.ind_a]->f_bh[ipol.ind_r[0]][ipol.ind_th[0]],
      tab->dat[ipol.ind_a]->f_bh[ipol.ind_r[0]][ipol.ind_th[0] + 1],
      tab->dat[ipol.ind_a]->f_bh[ipol.ind_r[0] + 1][ipol.ind_th[0]],
      tab->dat[ipol.ind_a]->f_bh[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1],
      tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_r[1]][ipol.ind_th[1]],
      tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_r[1]][ipol.ind_th[1] + 1],
      tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_r[1] + 1][ipol.ind_th[1]],
      tab->dat[ipol.ind_a + 1]->f_bh[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1]);

  fractions->refl_frac = fractions->f_ad / fractions->f_inf;

  if (fractions->f_inf > 0.5) { // photons are not allowed to cross the disk plane - should never be the case
    fractions->f_inf = 0.5;
    printf("Rellp_Extended: escaping photon fraction more than 0.5 (happens only for r_src > 40rg). \n");
    printf("Replacing with 0.5. This will be fixed in the new extended table. \n");
  }

  // these are equal for now. f_inf_rest is used to calculate primary spectrum
  fractions->f_inf_rest = fractions->f_inf; // true as long as beta=0

  return fractions;
}


double calc_primary_lensing(const double incl, const extendedSourceTable* prim_tab, const InterpFactors ipol) {
  int ind_incl;
  double ifac_incl;
  get_ipol_factor_double(incl, prim_tab->dat[ipol.ind_a]->bin, prim_tab->n_bin,
                         &ind_incl, &ifac_incl);
  if (ifac_incl < 0.0 or ifac_incl > 1.0) {
    // should not happen in the current table
    printf(" *** Warning : inclination interpolation becomes extrapolation: ifac_i = %.1f \n", ifac_incl);
    printf("Better limit a range of values within table limits, 0.05 < cos(incl) < 1.0 \n");
  }

  double lens_low = get_interpolated_value_a_r_th(
      ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0]][ind_incl],
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_incl],
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_incl],
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_incl],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1]][ind_incl],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_incl],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_incl],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_incl]);
  double lens_high = get_interpolated_value_a_r_th(
      ipol.ifac_a, ipol.ifac_r, ipol.ifac_th,
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0]][ind_incl + 1],
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0]][ipol.ind_th[0] + 1][ind_incl + 1],
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0]][ind_incl + 1],
      prim_tab->dat[ipol.ind_a]->val[ipol.ind_r[0] + 1][ipol.ind_th[0] + 1][ind_incl + 1],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1]][ind_incl + 1],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1]][ipol.ind_th[1] + 1][ind_incl + 1],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1]][ind_incl + 1],
      prim_tab->dat[ipol.ind_a + 1]->val[ipol.ind_r[1] + 1][ipol.ind_th[1] + 1][ind_incl + 1]);
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
void calc_ring_emis(emisProfile* emisProf, const double a, const double r_src,
                    const double theta, const double gamma, extendedSourceTable* tab, int* status) {

  CHECK_STATUS_VOID(*status);

  InterpFactors ipol;
  get_ipol_factors_a_r_th(a, r_src, theta, tab, ipol);

  double* radial_grid = nullptr;
  emisProfile *emis_profile_table = interpol_extended_table(tab, ipol, gamma, radial_grid, status);
  rebin_emisprofile_on_radial_grid(emisProf, emis_profile_table, status);

  free(radial_grid); // radial_grid is _not_ owned by emissivity profile structure
  free_emisProfile(emis_profile_table);


  // TODO: take the beaming of the jet into account (see Dauser et al., 2013)
//  for (int ii=0; ii<tab->n_rad; ii++){
//    if (param->beta > 1e-6) {
      // emisProf->emis[ii] *= pow(doppler_factor(emisProf->del_emit[ii], param->beta), 2);
//    }
//  }
}

/**
 * @brief Divides flux by Kerr-to-flat area ratio and by Lorentz factor
 * TODO: cache these corrections in new table, so this function will be obsolete
 */
void apply_flux_correction(extendedSourceData* dat, const double* r, const double spin,
                           const int n_r, const int n_th, const int n_rad) {
  for (int i_r = 0; i_r < n_r; i_r++) {
    for (int j_th = 0; j_th < n_th; j_th++) {
      for (int k_rad = 0; k_rad < n_rad; k_rad++){
        dat->val[i_r][j_th][k_rad] /= (kerr_to_flat_area_ratio(r[k_rad], spin) *
            gamma_correction(r[k_rad], spin));
      }
    }
  }
}

/**
 * @brief Calculates photon fractions (esc, bh, disk), lensing and boost of the primary emission for ring geometry
 * @details
 * @param emisProf; emisProfile, pointer to the structure with emissivity profile and necessary quantities
 * @param param; relParam, pointer to the parameters structure
 * @param prim_tab; extendedSourceTable, pointer to the primary table
 */
void calc_photon_fate_ring(emisProfile* emisProf, const relParam* param, extendedSourceTable* prim_tab, int* status) {
  CHECK_STATUS_VOID(*status);

  InterpFactors prim_ipol;
  get_ipol_factors_a_r_th(param->a, param->r_src, param->theta, prim_tab, prim_ipol);
  emisProf->photon_fate_fractions = calc_refl_frac_ext(prim_tab, prim_ipol, status);
  // input incl can be both in radian and in grads. The table has cos incl, so convert:
  double cos_incl = (param->incl < 0.5 * M_PI) ? cos(param->incl) : cos(param->incl * M_PI / 180.0);

  emisProf->photon_fate_fractions->lensing = calc_primary_lensing(cos_incl, prim_tab, prim_ipol);
  emisProf->photon_fate_fractions->lensing_and_boost_factor = (emisProf->photon_fate_fractions->lensing *
      calc_energy_shift_gamma(prim_tab, prim_ipol, cos_incl, param->gamma));
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
void calc_emis_ring_source(emisProfile* emisProf, const relParam* param, int* status) {

  CHECK_STATUS_VOID(*status);

  extendedSourceTable* emissivity_table = get_es_table(status);
  extendedSourceTable* prim_table = get_prim_table(status);

  calc_ring_emis(emisProf, param->a, param->r_src, param->theta, param->gamma, emissivity_table, status);
  calc_photon_fate_ring(emisProf, param, prim_table, status);

  CHECK_STATUS_VOID(*status);
}

void calc_emis_slab_source(emisProfile* emisProf, const relParam* param, int* status) {

  CHECK_STATUS_VOID(*status);

  extendedSourceTable* emissivity_table = get_es_table(status);
  extendedSourceTable* prim_table = get_prim_table(status);
  InterpFactors ipol;
  InterpFactors prim_ipol;
  // TODO proper scheme for primary normalization, currently we assume just the values from the outermost ring:
  calc_photon_fate_ring(emisProf, param, prim_table, status);

  const int num_rings = 100;
  std::array<double, num_rings> slab_grid{};
  emisProfile* slab_emis_sum = nullptr;
  double* rad_grid = nullptr;

  double weight_i;
  double weight_sum = 0.0;

  for (int ii = 0; ii < num_rings; ii++) { // define grid of rings for slab
    slab_grid[ii] = param->x_in + (param->x - param->x_in) * ii / (num_rings - 1);
  }

  for (int ii = 0; ii < num_rings; ii++) {

    /// !!! FIX for slab: either new function for interpolation in h,x plane, or convert slab rings to r, theta
    get_ipol_factors_a_r_th(param->a, param->height, slab_grid[ii], emissivity_table, ipol);
    get_ipol_factors_a_r_th(param->a, param->height, slab_grid[ii], prim_table, prim_ipol);

    emisProfile *ring_emis_tmp = interpol_extended_table(emissivity_table, ipol,
                                                         param->gamma, rad_grid, status);

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
      weight_i = slab_grid[ii] * (slab_grid[ii + 1] - slab_grid[ii]) * ring_emis_tmp->photon_fate_fractions->refl_frac;
    }
    weight_sum += weight_i;
    for (int nn = 0; nn < ring_emis_tmp->nr; nn++) {
      slab_emis_sum->emis[nn] += ring_emis_tmp->emis[nn] * weight_i;
    }
    // free(ring_emis_tmp->re);
    free_emisProfile(ring_emis_tmp);
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
 * @brief Gets g^Gamma from ring source to observer from the table (at a given inclination)
 * @details It is almost duplicated with the function below, which does the same for the disk point
 * @return The g^Gamma value with arbitrary Gamma (e.g., Gamma = 1 if you want just energy shift)
 */
double get_energyboost_ring_source_obs(const relParam* param, const double gamma) {
  int status = EXIT_SUCCESS; // in case it is loaded from somewhere in the code
  extendedSourceTable* prim_table = get_prim_table(&status);
  InterpFactors prim_ipol;
  get_ipol_factors_a_r_th(param->a, param->r_src, param->theta, prim_table, prim_ipol);
  // input incl can be both in radian and in grads. The table has cos incl, so convert it here to cos
  // (I use the fact that relxill limits are 3 and 87 deg, so if input incl is < 1.57 (~pi/2) - it is radian)
  double cos_incl = (param->incl < 0.5 * M_PI) ? cos(param->incl) : cos(param->incl * M_PI / 180.0);
  double gshift = calc_energy_shift_gamma(prim_table, prim_ipol, cos_incl, gamma);
  return gshift;
}

/**
 * @brief Gets g^Gamma from ring source to disk for a given disk radius
 * @param param relxill parameters and values
 * @param radius Disk radius
 * @param gamma Photon index
 * @return The g^Gamma value
 */
double get_energyboost_ring_source_disk(const relParam* param, const double radius, const double gamma) {
  int status = EXIT_SUCCESS;
  extendedSourceTable* table = get_es_table(&status);
  InterpFactors ipol;
  get_ipol_factors_a_r_th(param->a, param->r_src, param->theta, table, ipol);
  double gshift = calc_energy_shift_gamma(table, ipol, radius, gamma);
  return gshift;
}


// constructors and destructors //
/**
 * @brief Sets the memory for new extended source table
 * @details
 * @param n_r; int, number of spherical radii in the table
 * @param n_th; int, number of polar angles in the table
 * @param n_a; int, number of spins in the table
 * @param n_bin; int, number of bins in the table (radius or inclination)
 * @param n_storage; int, number of the energy shift bins in the table
 * @return new_extended_source_table, new extended source table
 */
extendedSourceTable* new_extended_source_table(const int n_r, const int n_th,
                                               const int n_a, const int n_bin, const int n_storage,
                                               int* status) {

  auto *tab = (extendedSourceTable*) malloc(sizeof(extendedSourceTable));

  CHECK_MALLOC_RET_STATUS(tab, status, nullptr)

  tab->n_a = n_a;
  tab->n_r = n_r;
  tab->n_bin = n_bin;
  tab->n_storage = n_storage;
  tab->n_th = n_th;

  tab->a = nullptr;

  tab->dat = nullptr;

  tab->dat = (extendedSourceData**) malloc(sizeof(extendedSourceData*) * tab->n_a);

  CHECK_MALLOC_RET_STATUS(tab->dat, status, tab)

  for (int i_a = 0; i_a < tab->n_a; i_a++) {
    tab->dat[i_a] = nullptr;
  }
  return tab;
}


/**
 *
 * @brief Sets the memory for new extended source data
 * @details
 * @return extendedSourceData, new extended source data
 */
extendedSourceData* new_extended_source_data(int* status, const bool is_lensing_table) {
  int n_r = EXT_TABLE_NUM_RADII;
  int n_th = EXT_TABLE_NUM_ANGLES;
  int n_bin = EXT_TABLE_NUM_RAD_BINS;
  int n_shift = EXT_TABLE_NUM_ENERGY_BINS;
  if (is_lensing_table) {
    n_r = PRIM_TABLE_NUM_RADII;
    n_th = PRIM_TABLE_NUM_ANGLES;
    n_bin = PRIM_TABLE_NUM_INCLINATIONS;
    n_shift = PRIM_TABLE_NUM_ENERGY_BINS;
  }
  auto *dat = (extendedSourceData*) malloc(sizeof(extendedSourceData));
  CHECK_MALLOC_RET_STATUS(dat, status, nullptr)

  dat->r_sph = (double*) malloc(sizeof(double) * n_r);

  CHECK_MALLOC_RET_STATUS(dat->r_sph, status, nullptr)

  dat->theta = (double**) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->theta, status, nullptr)

  dat->bin = (double*) malloc(sizeof(double) * n_bin);
  CHECK_MALLOC_RET_STATUS(dat->bin, status, nullptr)

  dat->refl_frac = (double**) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, nullptr)

  dat->f_ad = (double**) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->f_ad, status, nullptr)

  dat->f_inf = (double**) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->f_inf, status, nullptr)

  dat->f_bh = (double**) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->f_bh, status, nullptr)

  dat->val = (double***) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->val, status, nullptr)

  dat->gshift = (double***) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->gshift, status, nullptr)

  dat->energy_data_storage = (double****) malloc(sizeof(double) * n_r);
  CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, nullptr)

  for (int i_r = 0; i_r < n_r; i_r++) {

    dat->refl_frac[i_r] = (double*) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, nullptr)

    dat->theta[i_r] = (double*) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->theta, status, nullptr)

    dat->f_ad[i_r] = (double*) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->f_ad, status, nullptr)

    dat->f_inf[i_r] = (double*) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->f_inf, status, nullptr)

    dat->f_bh[i_r] = (double*) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->f_bh, status, nullptr)

    dat->gshift[i_r] = (double**) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->gshift, status, nullptr)

    dat->val[i_r] = (double**) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->val, status, nullptr)

    dat->energy_data_storage[i_r] = (double***) malloc(sizeof(double) * n_th);
    CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, nullptr)

    for(int j_th = 0; j_th < n_th; j_th++) {
      dat->gshift[i_r][j_th] = (double*) malloc(sizeof(double) * n_bin);
      dat->val[i_r][j_th] = (double*) malloc(sizeof(double) * n_bin);
      dat->energy_data_storage[i_r][j_th] = (double**) malloc(sizeof(double) * n_bin);

      for(int k_e = 0; k_e < n_bin; k_e++) {
        dat->energy_data_storage[i_r][j_th][k_e] = (double*) malloc(sizeof(double) * n_shift);
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
void free_extended_source_table(extendedSourceTable* tab, const bool is_lensing_table) {
  if (tab != nullptr) {
    if (tab->dat != nullptr) {
      for (int i_a = 0; i_a < tab->n_a; i_a++) {
        if (tab->dat[i_a] != nullptr){
          free_extended_source_data(tab->dat[i_a], is_lensing_table);
          free(tab->dat[i_a]);
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
void free_extended_source_data(extendedSourceData* dat, const bool is_lensing_table) {
  int n_r = EXT_TABLE_NUM_RADII;
  int n_th = EXT_TABLE_NUM_ANGLES;
  int n_bin = EXT_TABLE_NUM_RAD_BINS;
  if (is_lensing_table) {
    n_r = PRIM_TABLE_NUM_RADII; // re-define in case of lensing table, dimensions are different
    n_th = PRIM_TABLE_NUM_ANGLES;
    n_bin = PRIM_TABLE_NUM_INCLINATIONS;
  }
  if (dat != nullptr) {
    for (int i_r = 0; i_r < n_r; i_r++) {
      if (dat->theta != nullptr) free(dat->theta[i_r]);
      if (dat->val != nullptr) free(dat->val[i_r]);
      if (dat->gshift != nullptr) free(dat->gshift[i_r]);
      if (dat->refl_frac != nullptr) free(dat->refl_frac[i_r]);
      if (dat->f_ad != nullptr) free(dat->f_ad[i_r]);
      if (dat->f_inf != nullptr) free(dat->f_inf[i_r]);
      if (dat->f_bh != nullptr) free(dat->f_bh[i_r]);
      if (dat->energy_data_storage != nullptr) free(dat->energy_data_storage[i_r]);

      for (int j_th = 0; j_th < n_th; j_th++) {
        if (dat->val[i_r] != nullptr) free(dat->val[i_r][j_th]);
        if (dat->gshift[i_r] != nullptr) free(dat->gshift[i_r][j_th]);

        for (int k_b = 0; k_b < n_bin; k_b++) {
          if (dat->energy_data_storage[i_r][j_th] != nullptr) free(dat->energy_data_storage[i_r][j_th][k_b]);
        }
      }
    }
  }
    free(dat->bin);
    free(dat->r_sph);
    free(dat->theta);
    free(dat->val);
    free(dat->gshift);
    free(dat->energy_data_storage);
    free(dat->refl_frac);
    free(dat->f_inf);
    free(dat->f_ad);
    free(dat->f_bh);

}
