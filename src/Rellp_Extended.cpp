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
  bool lensing_table = true;
  if (cached_extendedSource_table_prim == nullptr) {
    read_extendedSource_table(PRIM_TABLE_FILENAME, &cached_extendedSource_table_prim, status,
                              lensing_table);
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
  // return value interpolated between spins
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
                               double spin, bool lensing) {

  extendedSourceData *dat = new_extendedSourceData(status, lensing);
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
      if (lensing) { // extra column which is now only in the primary table
        fits_read_col(fptr, TDOUBLE, 8, irow, 1, 1, &doublenullptr,
                      &dat->f_bh[ii][jj],&anynullptr, status);
      }
      dat->refl_frac[ii][jj] = dat->f_ad[ii][jj]/dat->f_inf[ii][jj];

    }
  }
  if (!lensing) { // only for flux, not for lensing
    /// to be removed with new table, which does these correction by default
    apply_flux_correction(dat, dat->bin, spin, n_h, n_r, n_bin);
  }
  return dat;
}


/**
 * @brief Calculates reflection fraction from interpolated table values, and photon fractions
 * @details
 * @param param, const relParam, a pointer to the structure of parameters used, the parameters include:
 * a - spin, height - height, r_ring - ring radius
 * @param tab, extendedSourceTable, a pointer to the structure containing a table
 * @return lpReflFrac, a pointer to the structure containing interpolated reflection fraction and photon fractions
 */
static lpReflFrac *calc_refl_frac_ext(const relParam *param,
                                      extendedSourceTable *tab, int *status) {

  // Note that this function partially repeats the interpol_estable function. Maybe combine them?
  // the difference is mainly the return.

  //// I could make a function "get interpolation factors" and pass the ifac_a, h, r etc to the relevant functions

  lpReflFrac *fractions = new_lpReflFrac(status);
  CHECK_STATUS_RET(*status, fractions);

  int ind_a;
  double ifac_a;
  get_ipol_factor((float) param->a, tab->a, tab->n_a, &ind_a, &ifac_a);
  extendedSourceData *dat_ind_a[2] = {tab->dat[ind_a], tab->dat[ind_a+1]};

  int ind_h[2];
  double ifac_h[2];

  for (int ii = 0; ii < 2; ii++) {
    get_ipol_factor_double(param->height, dat_ind_a[ii]->height, tab->n_h, &(ind_h[ii]),
                           &(ifac_h[ii]));
  }

  int ind_r[2];
  double ifac_r[2];
  for (int ii = 0; ii < 2; ii++) {
    int index;
    // doesn't matter for the currently used table (but if r_ring are different for different heights - it will matter)
    if (ifac_h[ii] > 0.999999) {
      index = ind_h[ii] + 1;
    } else {
      index = ind_h[ii];
    }
    get_ipol_factor_double(param->r_ring, dat_ind_a[ii]->r_ring[index], tab->n_r,
                           &(ind_r[ii]), &(ifac_r[ii]));
  }

  fractions->f_ad = get_interpolated_value(ifac_a, ifac_h, ifac_r,
      dat_ind_a[0]->f_ad[ind_h[0]][ind_r[0]],dat_ind_a[0]->f_ad[ind_h[0]][ind_r[1]],
      dat_ind_a[0]->f_ad[ind_h[1]][ind_r[0]],dat_ind_a[0]->f_ad[ind_h[1]][ind_r[1]],
      dat_ind_a[1]->f_ad[ind_h[0]][ind_r[0]],dat_ind_a[1]->f_ad[ind_h[0]][ind_r[1]],
      dat_ind_a[1]->f_ad[ind_h[1]][ind_r[0]],dat_ind_a[1]->f_ad[ind_h[1]][ind_r[1]]);

  fractions->f_inf = get_interpolated_value(ifac_a, ifac_h, ifac_r,
      dat_ind_a[0]->f_inf[ind_h[0]][ind_r[0]],dat_ind_a[0]->f_inf[ind_h[0]][ind_r[1]],
      dat_ind_a[0]->f_inf[ind_h[1]][ind_r[0]],dat_ind_a[0]->f_inf[ind_h[1]][ind_r[1]],
      dat_ind_a[1]->f_inf[ind_h[0]][ind_r[0]],dat_ind_a[1]->f_inf[ind_h[0]][ind_r[1]],
      dat_ind_a[1]->f_inf[ind_h[1]][ind_r[0]],dat_ind_a[1]->f_inf[ind_h[1]][ind_r[1]]);
  // new table contains this value too
  fractions->f_bh = get_interpolated_value(ifac_a, ifac_h, ifac_r,
      dat_ind_a[0]->f_bh[ind_h[0]][ind_r[0]],dat_ind_a[0]->f_bh[ind_h[0]][ind_r[1]],
      dat_ind_a[0]->f_bh[ind_h[1]][ind_r[0]],dat_ind_a[0]->f_bh[ind_h[1]][ind_r[1]],
      dat_ind_a[1]->f_bh[ind_h[0]][ind_r[0]],dat_ind_a[1]->f_bh[ind_h[0]][ind_r[1]],
      dat_ind_a[1]->f_bh[ind_h[1]][ind_r[0]],dat_ind_a[1]->f_bh[ind_h[1]][ind_r[1]]);

  // before there was one more interpolation of refl_frac between values written in function loadSingle source ... above
  // but I don't see any obstacles against just calculating it here. Need to check if it is calculated the same
  fractions->refl_frac = fractions->f_ad / fractions->f_inf;
  // photons are not allowed to cross the disk plane
  if (fractions->f_inf > 0.5) { // should never be the case - TODO in new table
    fractions->f_inf = 0.5;
    printf("Rellp_Extended: inf fraction more than 0.5. \n");
    printf("Replacing with 0.5. This will be fixed in the new extended table. \n");
  }
  /// these are equal for now. f_inf_rest is used to calculate primary spectrum
  /// (function flux_boost_source_to_observer in PrimarySource.h)
  /// later maybe would be good to distinguish them correctly - AN
  fractions->f_inf_rest = fractions->f_inf; /// anyway this is true for non-moving sources (beta=0)
  // str->f_inf_rest = 0.5 * (1.0 + cos(del_emit_ad_max));

  ////str->f_bh = 1.0 - (str->f_ad + str->f_inf); // this is not true if we re-define f_inf to be <= 0.5.
  // But it is important only for very large heights.
  // TODO: This definition does not take into account a possible variability of the disk edges. Is there any way to
  // overcome this limitation? Probably only approximate way, because with raytracing it is not easy to vary r_in/r_out

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

  extendedSourceTable* emissivity_table = get_es_table(status); // reading table everytime? No, caching
  // new part - loading primary normalization, i.e., lensing * g_{so}^Gamma, to primary normalization
  extendedSourceTable* prim_table = get_prim_table(status);
  calc_extended_emis(emis, param, emissivity_table, prim_table, status); // PRIM , prim_table

}


/**
 * @brief reads extended source table
 * @details
 * @param filename; const char, pointer to the name of the file (now declared in Rellp_Extended.h)
 * @param inp_tab; extendedSourceTable, two-fold pointer to the input table data
 */
void read_extendedSource_table(const char *filename, extendedSourceTable **inp_tab, int *status, bool lensing_table_flag) {

  extendedSourceTable *tab = (*inp_tab);
  fitsfile *fptr = nullptr;

  char fullfilename[999];

  do {
    if (tab != nullptr) {
      RELXILL_ERROR("extended source table already loaded\n", status);
      break;
    }
    if (lensing_table_flag) {
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
          tab->a[ii], lensing_table_flag);
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
    free_extendedSourceTable(tab, lensing_table_flag);
  }

  if(fptr != nullptr) {fits_close_file(fptr, status);}

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
////////////////////////////////////////////////////////
// check this func
  emisProfile *emis_profile_table = interpol_estable(
      param->a, param->height, param->r_ring, param->gamma, tab, status);
  // what about rebining of lensing?
   rebin_emisprofile_on_radial_grid(emisProf, emis_profile_table, status);

  free(emis_profile_table->re); // is not freed by free_emisProfile (Rellp.cpp)
  free_emisProfile(emis_profile_table);

// check this func
  //PRIM emisProf->photon_fate_fractions = calc_refl_frac_ext(param, tab_prim, status);
  emisProf->photon_fate_fractions = calc_refl_frac_ext(param, prim_tab, status);

  // emisProf->normFactorPrimSpec = calc_primary_lensing(param, tab_prim, status);
  //// emisProf->normFactorPrimSpec *= calc_energy_shift_gamma(param, tab_prim, status);


  // take the beaming of the jet into account (see Dauser et al., 2013)
  // Not done for extended source, todo later?
//  for (int ii=0; ii<tab->n_rad; ii++){
//    if (param->beta > 1e-6) {
//      emisProf->emis[ii] *= pow(doppler_factor(emisProf->del_emit[ii], param->beta), 2);
//    }
//  }
 // AN: not used as far as I see from other routines?
 // emisProf->normFactorPrimSpec = calc_norm_factor_primary_spectrum_ext(param->height, param->a, param->gamma,
 //         emisProf->photon_fate_fractions->f_inf);
}

/**
 * @brief determine the gshift^Gamma for table grid values. Required for further interpolation
 */
void calc_energy_shift_gamma_tabulated(extendedSourceData &dat_ind_a, const int ind_h[2], const int ind_r[2],
                             double ****energy_shift_bins, int num_radial_bins, int num_shift_bins, double gamma){
  for (int i = 0; i < 2; i++) { // vertical bins, ind_h
    for (int j = 0; j < 2; j++) { // radial extent bins, ind_r
      for (int k = 0; k < num_radial_bins; k++) { // disk radial bins
        double sum = 0.0; // sum is zeroed every time
        for (int l = 0; l < num_shift_bins; l++) { // energy shift bins
          sum += pow(energy_shift_bins[ind_h[i]][ind_r[j]][k][l], gamma);
          assert(sum >= 0.0);
        }
        // replace with -> if will do dynamical allocation
        dat_ind_a.gshift[ind_h[i]][ind_r[j]][k] = sum / num_shift_bins;
      }
    }
  }
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
emisProfile *interpol_estable(double a, double height, double r_ring,
                              double gamma, extendedSourceTable *tab, int *status) {
  CHECK_STATUS_RET(*status, nullptr);

  int ind_a;
  double ifac_a;
  get_ipol_factor((float) a, tab->a, tab->n_a, &ind_a, &ifac_a);

  // assert(ifac_a >= 0.0 && ifac_a <= 1.0);
  if (ifac_a < 0.0 or ifac_a > 1.0){
    // should not happen in the current table
    printf(" *** Warning : spin interpolation becomes extrapolation: ifac_a = %.1f \n", ifac_a);
    printf("Better limit a range of values within table limits \n");
  }
  extendedSourceData *dat_spin[2] = {tab->dat[ind_a], tab->dat[ind_a+1]};

  // Need to interpolate radial grid of table because it is different for different spins (risco(a))
  auto* radial_grid = (double*) malloc(sizeof(double)*tab->n_bin);
  CHECK_MALLOC_RET_STATUS(radial_grid,status,nullptr)
  for (int ii = 0; ii < tab->n_bin; ii++) {
    radial_grid[ii] = interp_lin_1d(ifac_a, dat_spin[0]->bin[ii], dat_spin[1]->bin[ii]);
  }

  emisProfile* emis_profile_table = new_emisProfile(radial_grid, tab->n_bin, status);

  int ind_h[2];
  double ifac_h[2];

  for (int ii = 0; ii < 2; ii++) {
    // interpolation factors for height
    get_ipol_factor_double(height, dat_spin[ii]->height, tab->n_h, &(ind_h[ii]),
                           &(ifac_h[ii]));

    if (ifac_h[ii] < 0.0 or ifac_h[ii] > 1.0){
      printf(" *** Warning : height interpolation becomes extrapolation: ifac_h = %.1f \n", ifac_h[ii]);
      printf("Better limit a range of values within table limits \n");
    }
  }

  int ind_r[2];
  double ifac_r[2];
  for (int ii = 0; ii < 2; ii++) {
    int radial_index; // index of the range of ring values corresponding to a given height
    // (in the currently used table it does not matter, as ring values in table are equal for all heights)
    // if h almost equals the upper height, we also shift the index of the ring value
    if (ifac_h[ii] > 0.999999) {
      radial_index = ind_h[ii] + 1;
    } else {
      radial_index = ind_h[ii];
    }
    // interpolation factors for extent
    get_ipol_factor_double(r_ring, dat_spin[ii]->r_ring[radial_index], tab->n_r, &(ind_r[ii]),
                           &(ifac_r[ii]));
    if (ifac_r[ii] < 0.0 or ifac_r[ii] > 1.0){
      printf(" *** Warning : dist interpolation becomes extrapolation: ifac_r = %.1f \n", ifac_r[ii]);
      printf("Better limit a range of values within table limits \n");
    }
  }

  calc_energy_shift_gamma_tabulated(*dat_spin[0], ind_h, ind_r, tab->dat[ind_a]->energy_data_storage,
                          tab->n_bin, tab->n_storage, gamma);
  calc_energy_shift_gamma_tabulated(*dat_spin[1], ind_h, ind_r, tab->dat[ind_a + 1]->energy_data_storage,
                          tab->n_bin, tab->n_storage, gamma);

  for (int ii = 0; ii < tab->n_bin; ii++) {
    // Interpolation of the flux
    // looks like meh, but gives some generalization potential (maybe refactor to smth better reading further?)
    emis_profile_table->emis[ii] = get_interpolated_value(ifac_a, ifac_h, ifac_r,
        dat_spin[0]->val[ind_h[0]][ind_r[0]][ii], dat_spin[0]->val[ind_h[0]][ind_r[1]][ii],
        dat_spin[0]->val[ind_h[1]][ind_r[0]][ii], dat_spin[0]->val[ind_h[1]][ind_r[1]][ii],
            dat_spin[1]->val[ind_h[0]][ind_r[0]][ii], dat_spin[1]->val[ind_h[0]][ind_r[1]][ii],
            dat_spin[1]->val[ind_h[1]][ind_r[0]][ii], dat_spin[1]->val[ind_h[1]][ind_r[1]][ii]);

    /// would be better if it was a separate function taking current rel_params as input. Need it for normalization
    double g_gamma = get_interpolated_value(ifac_a, ifac_h, ifac_r,
        dat_spin[0]->gshift[ind_h[0]][ind_r[0]][ii], dat_spin[0]->gshift[ind_h[0]][ind_r[1]][ii],
        dat_spin[0]->gshift[ind_h[1]][ind_r[0]][ii], dat_spin[0]->gshift[ind_h[1]][ind_r[1]][ii],
            dat_spin[1]->gshift[ind_h[0]][ind_r[0]][ii], dat_spin[1]->gshift[ind_h[0]][ind_r[1]][ii],
            dat_spin[1]->gshift[ind_h[1]][ind_r[0]][ii], dat_spin[1]->gshift[ind_h[1]][ind_r[1]][ii]);
    // printf("check emis and gGamma, %.10f, %.10f\n", emis_profile_table->emis[ii], g_gamma);
    emis_profile_table->emis[ii] *= g_gamma; // this line accounts for g^Gamma

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
void free_extendedSourceData(extendedSourceData *dat, bool lensing) {
  int n_h = EXT_TABLE_NHEIGHT;
  int n_x = EXT_TABLE_NDISTANCE;
  int n_r = EXT_TABLE_NRAD;
  if (lensing) {
    n_h = PRIM_TABLE_NHEIGHT; // re-define in case of lensing table, dimensions can be different
    n_x = PRIM_TABLE_NDISTANCE;
    n_r = PRIM_TABLE_NINCL;
  }
  if (dat != nullptr) {
    for (int ii =0; ii < n_h; ii++) {
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

}
