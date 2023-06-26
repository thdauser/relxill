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
 * @brief This module creates emissivity profile ready for usage in relxill taking it from the input table
 * @details This module reads values from the table produced with corona-lean by Stefan Licklederer.
 * The table is read by the columns, the columns are interpolated for further usage.
 */


#include "cfitsio.h"
// #include "stdio.h"
#include "Relbase.h"
#include "Rellp_Extended.h"
#include "Rellp.h"
#include "Relphysics.h"

extern "C"{
#include "relutility.h"
#include "writeOutfiles.h"
}

extendedSourceTable *cached_extendedSource_table = nullptr;

/**
 * @brief Gets extended source table with a given name
 * @details
 * @param status pointer, status of the code for debugging
 * @return extendedSourceTable, cached extended source table
 */
extendedSourceTable* get_es_table(int* status){
  CHECK_STATUS_RET(*status,nullptr);

  if (cached_extendedSource_table == nullptr) {
    read_extendedSource_table(ESTABLE_FILENAME, &cached_extendedSource_table, status);
    CHECK_STATUS_RET(*status, nullptr);
  }
  return cached_extendedSource_table;
}

/**
 * @brief Gets spin values from the table for each row
 * @details Uses cfitsio to read preloaded fits file
 * @param nrows integer, number of rows in the table
 * @param val two-fold pointer, to the value to be read (or what?)
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
 * @brief Gets data for a single source from the table and then makes a relativistic correction of the flux
 * @details Need to set the variables for fits_read_col instead of using numbers with unclear origin
 * Should the relativistic correction be made by a separate function better?
 * @param fptr pointer, to the fits file location
 * @param n_h integer, number of heights simulated
 * @param n_distance integer, number of distances simulated
 * @param n_rad integer, number of (radial bins??)
 * @param n_storage integer, number of (what a storage??). Used only to read the fifth column
 * @param spin double, spin value of the source
 * @return extendedSourceData, data from the table
 */
static extendedSourceData *
load_single_extendedSourceData(fitsfile *fptr, int n_h, int n_distance,
                               int n_rad, int n_storage, int *status,
                               double spin) {

  extendedSourceData *dat = new_extendedSourceData(status);
  CHECK_MALLOC_RET_STATUS(dat, status, nullptr)


  int anynullptr = 0;
  double doublenullptr = 0.0;

  /*  read the columns */
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, n_rad, &doublenullptr, dat->radius,
                &anynullptr, status);

  for (int ii = 0; ii < n_h; ii++) {
    fits_read_col(fptr, TDOUBLE, 2, 1 + ii * n_h, 1, 1, &doublenullptr, &(dat->height[ii]),
                  &anynullptr, status);
    fits_read_col(fptr, TDOUBLE, 3, 1 + ii * n_h, 1, n_distance, &doublenullptr, dat->distance[ii],
                  &anynullptr, status);
  }

  for (int ii = 0; ii < n_h; ii++) {
    for (int jj = 0; jj < n_distance; jj++) {

      int irow = (ii * n_distance) + jj + 1;

      fits_read_col(fptr, TDOUBLE, 4, irow, 1, n_rad, &doublenullptr, dat->flux[ii][jj],
                    &anynullptr, status);

      for (int kk = 0; kk < n_rad; kk++){
        fits_read_col(fptr, TDOUBLE, 5, irow, kk*n_storage+1, n_storage, &doublenullptr, dat->energy_data_storage[ii][jj][kk],
                      &anynullptr, status);
      }

      // n_rad to 1 because it is not an array just a single value
      fits_read_col(fptr, TDOUBLE, 6, irow, 1, 1, &doublenullptr, &dat->f_ad[ii][jj],
                    &anynullptr, status);

      fits_read_col(fptr, TDOUBLE, 7, irow, 1, 1, &doublenullptr, &dat->f_inf[ii][jj],
                    &anynullptr, status);

      dat->refl_frac[ii][jj] = dat->f_ad[ii][jj]/dat->f_inf[ii][jj];

    }
  }

  ////  add correction to flux
  double *r = dat->radius;
  double result;
  for (int i=0; i<n_h; i++) {
    for (int j=0; j<n_distance; j++) {
      for (int k=0; k<n_rad; k++){
        result = area_kerr(r[k],spin)*gamma_correction(r[k],spin);
        dat->flux[i][j][k] *= area_flat(r[k])/result;
      }
    }
  }

  return dat;
}


/**
 * @brief Modifies (?) a table tab making an interpolation of its content to the grid used in further calculations
 * @details
 * @param param, const relParam, a pointer to the structure of parameters used, the parameters include:
 * a - spin, height - height, distance - distance
 * @param tab, extendedSourceTable, a pointer to the structure containing a table
 * @return lpReflFrac, a pointer to the structure containing some doubles with parameters of bh and system (?)
 * probably the interpolated ones
 */
static lpReflFrac *calc_refl_frac_ext(const relParam *param,
                                      extendedSourceTable *tab, int *status) {

  // in case there is no relativity information, the refl_frac is 1

  lpReflFrac *str = new_lpReflFrac(status);
  CHECK_STATUS_RET(*status, str);


  // BEGINNING OF REFL_FRAC INTERPOLATION
  int ind_a;
  double ifac_a;
  get_ipol_factor((float) param->a, tab->a, tab->n_a, &ind_a, &ifac_a);
  extendedSourceData *dat_ind_a[2] = {tab->dat[ind_a], tab->dat[ind_a+1]};

  int ind_h[2];
  double ifac_h[2];

  for (int ii = 0; ii < 2; ii++) {
    get_ipol_factor_double(param->height, dat_ind_a[ii]->height, tab->n_h, &(ind_h[ii]), &(ifac_h[ii]));
  }

  int ind_d[2];
  double ifac_d[2];
  for (int ii = 0; ii < 2; ii++) {
    int index;
    // ad hoc solution to solve issue with false choosing of distance value
    if (ifac_h[ii] > 0.999999) {
      index = ind_h[ii] + 1;
    } else {
      index = ind_h[ii];
    }
    get_ipol_factor_double(param->distance, dat_ind_a[ii]->distance[index], tab->n_offaxis, &(ind_d[ii]), &(ifac_d[ii]));
  }


  // interpolate f_ad value
  double rlo =
      interp_lin_2d(ifac_d[0], // does h-d interpolation for spin[ind_a]
                    ifac_h[0], dat_ind_a[0]->f_ad[ind_h[0]][ind_d[0]],
                    dat_ind_a[0]->f_ad[ind_h[0]][ind_d[0] + 1],
                    dat_ind_a[0]->f_ad[ind_h[0] + 1][ind_d[0]],
                    dat_ind_a[0]->f_ad[ind_h[0] + 1][ind_d[0] + 1]);
  double rhi = interp_lin_2d(
      ifac_d[1], // does height-distance interpolation for spin[ind_a+1]
      ifac_h[1], dat_ind_a[1]->f_ad[ind_h[1]][ind_d[1]],
      dat_ind_a[1]->f_ad[ind_h[1]][ind_d[1] + 1],
      dat_ind_a[1]->f_ad[ind_h[1] + 1][ind_d[1]],
      dat_ind_a[1]->f_ad[ind_h[1] + 1][ind_d[1] + 1]);
  str->f_ad = interp_lin_1d(ifac_a, rlo, rhi);

  // interpolate f_inf value
  double rlo_1 =
      interp_lin_2d(ifac_d[0], // does h-d interpolation for spin[ind_a]
                    ifac_h[0], dat_ind_a[0]->f_inf[ind_h[0]][ind_d[0]],
                    dat_ind_a[0]->f_inf[ind_h[0]][ind_d[0] + 1],
                    dat_ind_a[0]->f_inf[ind_h[0] + 1][ind_d[0]],
                    dat_ind_a[0]->f_inf[ind_h[0] + 1][ind_d[0] + 1]);
  double rhi_1 = interp_lin_2d(
      ifac_d[1], // does height-distance interpolation for spin[ind_a+1]
      ifac_h[1], dat_ind_a[1]->f_inf[ind_h[1]][ind_d[1]],
      dat_ind_a[1]->f_inf[ind_h[1]][ind_d[1] + 1],
      dat_ind_a[1]->f_inf[ind_h[1] + 1][ind_d[1]],
      dat_ind_a[1]->f_inf[ind_h[1] + 1][ind_d[1] + 1]);
  str->f_inf = interp_lin_1d(ifac_a, rlo_1, rhi_1);

  // interpolate refl_frac value
  double rlo_2 =
        interp_lin_2d(ifac_d[0], // does h-d interpolation for spin[ind_a]
                      ifac_h[0], dat_ind_a[0]->refl_frac[ind_h[0]][ind_d[0]],
                      dat_ind_a[0]->refl_frac[ind_h[0]][ind_d[0] + 1],
                      dat_ind_a[0]->refl_frac[ind_h[0] + 1][ind_d[0]],
                      dat_ind_a[0]->refl_frac[ind_h[0] + 1][ind_d[0] + 1]);
  double rhi_2 = interp_lin_2d(
        ifac_d[1], // does height-distance interpolation for spin[ind_a+1]
        ifac_h[1], dat_ind_a[1]->refl_frac[ind_h[1]][ind_d[1]],
        dat_ind_a[1]->refl_frac[ind_h[1]][ind_d[1] + 1],
        dat_ind_a[1]->refl_frac[ind_h[1] + 1][ind_d[1]],
        dat_ind_a[1]->refl_frac[ind_h[1] + 1][ind_d[1] + 1]);
  str->refl_frac = interp_lin_1d(ifac_a, rlo_2, rhi_2);
  // photons are not allowed to cross the disk plane
  if (str->f_inf > 0.5) {
    str->f_inf = 0.5;
  }

  str->f_bh = 1 - (str->f_ad + str->f_inf);

  return str;
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

  extendedSourceTable* estable = get_es_table(status);
  calc_extended_emis(emis, param, estable, status);

}


/**
 * @brief reads extended source table
 * @details
 * @param filename; const char, pointer to the name of the file (now declared in Rellp_Extended.h)
 * @param inp_tab; extendedSourceTable, two-fold pointer to the input table data
 */
void read_extendedSource_table(const char *filename, extendedSourceTable **inp_tab, int *status) {

  extendedSourceTable *tab = (*inp_tab);
  fitsfile *fptr = nullptr;

  char fullfilename[999]; // not used now but it would be good to add this option - writing local file name instead of full name

  do {
    if (tab != nullptr) {
      RELXILL_ERROR("extended source table already loaded\n", status);
      break;
    }

    tab = new_extendedSourceTable(ESTABLE_NHEIGHT, ESTABLE_NRAD,
                                  ESTABLE_NDISTANCE, ESTABLE_NA, ESTABLE_NENERGYBINS, status);
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
    printf("************************\n");
    printf("filename : %s\n", fullfilename);


    get_extendedsourcetable_spin(tab->n_a, &(tab->a) , fptr, status);
    CHECK_RELXILL_ERROR("reading of spin axis failed\n", status);


    // now load the full table (need to go through all extensions
    int ii;
    int rownum;
    for (ii = 0; ii < tab->n_a; ii++) {

      rownum = ii + 2; // number of the row we read in the fits table
      fits_movabs_hdu(fptr, rownum, nullptr, status);

      assert(tab->dat[ii] == nullptr);

      tab->dat[ii] = load_single_extendedSourceData(
          fptr, tab->n_h, tab->n_offaxis, tab->n_rad, tab->n_storage, status,
          tab->a[ii]);
      if (*status != EXIT_SUCCESS) {
        RELXILL_ERROR("failed to load data from the extended source table into memory\n", status);

        break;
      }
    }

  } while (0);

  if (*status == EXIT_SUCCESS) {
    // assigne the value
    (*inp_tab) = tab;
  } else {
    free_extendedSourceTable(tab);
  }

  if(fptr != nullptr) {fits_close_file(fptr, status);}

}


/////////////////////////////////////////////////////////
/**
 * @brief Calculates emissivity
 * @details
 * @param emisProf; emisProfile, pointer to the structure with emissivity profile
 * @param param; relParam, pointer to the parameters structure, as declared above
 * @param tab; extendedSourceTable, pointer to the table structure, as declared above
 */
void calc_extended_emis(emisProfile *emisProf, const relParam *param,
                        extendedSourceTable *tab, int *status) {

  CHECK_STATUS_VOID(*status);

  emisProfile* emis_profile_table = interpol_estable(
      param->a, param->height, param->d_ring, param->gamma, ESTABLE_NENERGYBINS, tab, status);

  rebin_emisprofile_on_radial_grid(emisProf, emis_profile_table, status);

  free(emis_profile_table->re); // is not freed by free_emisProfile
  free_emisProfile(emis_profile_table);

  emisProf->photon_fate_fractions = calc_refl_frac_ext(param, tab, status);

  // take the beaming of the jet into account (see Dauser et al., 2013)
  // Maybe not necessary for extended source
//  for (int ii=0; ii<tab->n_rad; ii++){
//    if (param->beta > 1e-6) {
//      emisProf->emis[ii] *= pow(doppler_factor(emisProf->del_emit[ii], param->beta), 2);
//    }
//  }

  emisProf->normFactorPrimSpec =
      calc_norm_factor_primary_spectrum_ext(param->height, param->a, param->gamma, emisProf->photon_fate_fractions->f_inf);
}


/**
 * @brief Interpolation of the table (not clear yet)
 * @details
 * @param a; double, spin
 * @param height; double, height
 * @param distance; double, distance
 * @param gamma; double, photon index
 * @param n_gshift_storage; int, number of energy shift bins
 * @param tab; extendedSourceTable, pointer to the table as declared above
 * @return emis_profile_table; emisProfile, pointer to the structure with emissivity profile
 */
emisProfile *interpol_estable(double a, double height, double distance,
                              double gamma, int n_gshift_storage,
                              extendedSourceTable *tab, int *status) {
  CHECK_STATUS_RET(*status, nullptr);

  int ind_a;
  double ifac_a;
  get_ipol_factor((float) a, tab->a, tab->n_a, &ind_a, &ifac_a);

  extendedSourceData *dat_ind_a[2] = {tab->dat[ind_a], tab->dat[ind_a+1]};

  double* jet_rad = (double*) malloc(sizeof(double)*tab->n_rad);
  CHECK_MALLOC_RET_STATUS(jet_rad,status,nullptr)
  for (int ii = 0; ii < tab->n_rad; ii++) {
    jet_rad[ii] = interp_lin_1d(ifac_a, dat_ind_a[0]->radius[ii], dat_ind_a[1]->radius[ii]);
  }

  emisProfile* emis_profile_table = new_emisProfile(jet_rad, tab->n_rad, status);

  int ind_h[2];
  double ifac_h[2];

  for (int ii = 0; ii < 2; ii++) {
    get_ipol_factor_double(height, dat_ind_a[ii]->height, tab->n_h, &(ind_h[ii]), &(ifac_h[ii]));
  }

  int ind_d[2];
  double ifac_d[2];
  for (int ii = 0; ii < 2; ii++) {
    int index;
    // ad hoc solution to solve issue with false choosing of distance value
    if (ifac_h[ii] > 0.999999) {
      index = ind_h[ii] + 1;
    } else {
      index = ind_h[ii];
    }
    get_ipol_factor_double(distance, dat_ind_a[ii]->distance[index], tab->n_offaxis, &(ind_d[ii]), &(ifac_d[ii]));
  }

  // determines the gshift for the two spin values required for
  // interpolation in the next function
  for(int i=0; i < tab->n_h; i++){
    for (int j=0; j< tab-> n_offaxis ; j++){
      for(int k=0; k < tab->n_rad ; k++){
        double sum1 = 0.0;
        double sum2 = 0.0;
        for(int l=0; l < tab->n_storage ; l++){
          sum1 += pow(tab->dat[ind_a]->energy_data_storage[i][j][k][l] / n_gshift_storage, gamma);
          sum2 += pow(tab->dat[ind_a+1]->energy_data_storage[i][j][k][l] / n_gshift_storage, gamma);
        }
        tab->dat[ind_a]->gshift[i][j][k] = sum1;
        tab->dat[ind_a+1]->gshift[i][j][k] = sum2;
      }
    }
  }

  interpol_emisprofile_spin_height_distance(emis_profile_table,
                                            ifac_a,
                                            ifac_h,
                                            ifac_d,
                                            ind_h,
                                            ind_d,
                                            dat_ind_a,
                                            tab->n_rad);


  return emis_profile_table;
}


/**
 *
 * @brief Interpolates values in height-distance plane
 * @details (if I understand correctly)
 * @param emis_profile_table; emisProfile, pointer to the structure containing emissivity profile table
 * @param ifac_a; double, interpolated? factor? for spin
 * @param ifac_h; const double, pointer to interpolated? factor? for height
 * @param ifac_d; const double, pointer to interpolated? factor? for distance
 * @param ind_h; const int, pointer to interpolated? index? for height
 * @param ind_d; const int, pointer to interpolated? index? for distance
 * @param dat_a; *const extendedSourceData (?), pointer to the structure with data
 * @param nr; const int, number of ?radial bins?
 */
void interpol_emisprofile_spin_height_distance(emisProfile *emis_profile_table,
                                               double ifac_a,
                                               const double *ifac_h,
                                               const double *ifac_d,
                                               const int *ind_h,
                                               const int *ind_d,
                                               extendedSourceData *const *dat_a,
                                               const int nr) {

  for (int ii = 0; ii < nr; ii++) {
    // Interpolation for the flux
    double rlo = interp_lin_2d(ifac_d[0], // does h-d interpolation for spin[ind_a]
                               ifac_h[0],
                               dat_a[0]->flux[ind_h[0]][ind_d[0]][ii],
                               dat_a[0]->flux[ind_h[0]][ind_d[0]+1][ii],
                               dat_a[0]->flux[ind_h[0]+1][ind_d[0]][ii],
                               dat_a[0]->flux[ind_h[0]+1][ind_d[0]+1][ii]);
    double rhi = interp_lin_2d(ifac_d[1], // does height-distance interpolation for spin[ind_a+1]
                               ifac_h[1],
                               dat_a[1]->flux[ind_h[1]][ind_d[1]][ii],
                               dat_a[1]->flux[ind_h[1]][ind_d[1]+1][ii],
                               dat_a[1]->flux[ind_h[1]+1][ind_d[1]][ii],
                               dat_a[1]->flux[ind_h[1]+1][ind_d[1]+1][ii]);
    emis_profile_table->emis[ii] = interp_lin_1d(ifac_a, rlo, rhi);

    double flo = interp_lin_2d(ifac_d[0], // does h-d interpolation for spin[ind_a]
                               ifac_h[0],
                               dat_a[0]->gshift[ind_h[0]][ind_d[0]][ii],
                               dat_a[0]->gshift[ind_h[0]][ind_d[0]+1][ii],
                               dat_a[0]->gshift[ind_h[0]+1][ind_d[0]][ii],
                               dat_a[0]->gshift[ind_h[0]+1][ind_d[0]+1][ii]);
    double fhi = interp_lin_2d(ifac_d[1], // does height-distance interpolation for spin[ind_a+1]
                               ifac_h[1],
                               dat_a[1]->gshift[ind_h[1]][ind_d[1]][ii],
                               dat_a[1]->gshift[ind_h[1]][ind_d[1]+1][ii],
                               dat_a[1]->gshift[ind_h[1]+1][ind_d[1]][ii],
                               dat_a[1]->gshift[ind_h[1]+1][ind_d[1]+1][ii]);

    // HOWEVER THE ENERGY SHIFT IS ALREADY TAKEN TO THE POWER OF GAMMA
    // NO NEED TO RAISE TO THE POWER HERE
    emis_profile_table->emis[ii] *= interp_lin_1d(ifac_a, flo, fhi);

    emis_profile_table->del_emit[ii] = 0; // TODO: Just as a placeholder to initialize value. Will remove once calculation is implemented.
    emis_profile_table->del_inc[ii] = 30 ;
  }
}

// constructors and desctructors //
/**
 *
 * @brief Sets the memory for new extended source table
 * @details
 * @param n_h; int, number of heights in the table
 * @param n_rad; int, number of bins in the table
 * @param n_offaxis; int, number of off-axis sources (??) in the table
 * @param n_a; int, number of spins in the table
 * @param n_storage; int, number of storage (??) in the table
 * @return new_extendedSourceTable, new extended source table
 */
extendedSourceTable *new_extendedSourceTable(int n_h, int n_rad, int n_offaxis,
                                             int n_a, int n_storage,
                                             int *status) {
  extendedSourceTable *tab = (extendedSourceTable *) malloc(sizeof(extendedSourceTable));

  CHECK_MALLOC_RET_STATUS(tab, status, nullptr)

  tab->n_a = n_a;
  tab->n_h = n_h;
  tab->n_rad = n_rad;
  tab->n_storage = n_storage;
  tab->n_offaxis = n_offaxis;

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
extendedSourceData *new_extendedSourceData(int *status) {

  extendedSourceData *dat = (extendedSourceData *) malloc(sizeof(extendedSourceData));
  CHECK_MALLOC_RET_STATUS(dat, status, nullptr)

  dat->height = (double *) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->height, status, nullptr)

  dat->distance = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->distance, status, nullptr)

  dat->radius = (double *) malloc(sizeof(double) * ESTABLE_NRAD);
  CHECK_MALLOC_RET_STATUS(dat->radius, status, nullptr)

  dat->refl_frac = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, nullptr)

  dat->f_ad = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->f_ad, status, nullptr)

  dat->f_inf = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->f_inf, status, nullptr)

  dat->flux = (double ***) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->flux, status, nullptr)

  dat->gshift = (double ***) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->flux, status, nullptr)

  dat->energy_data_storage = (double ****) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, nullptr)

  for (int i=0; i<ESTABLE_NHEIGHT; i++){

    dat->refl_frac[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, nullptr)

    dat->distance[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->distance, status, nullptr)

    dat->f_ad[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->f_ad, status, nullptr)

    dat->f_inf[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->f_inf, status, nullptr)

    dat->gshift[i] = (double **) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->gshift, status, nullptr)

    dat->flux[i] = (double **) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->flux, status, nullptr)

    dat->energy_data_storage[i] = (double ***) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, nullptr)

    for(int j=0; j<ESTABLE_NDISTANCE; j++){
      dat->gshift[i][j] = (double *) malloc(sizeof(double) * ESTABLE_NRAD);
      dat->flux[i][j] = (double *) malloc(sizeof(double) * ESTABLE_NRAD);
      dat->energy_data_storage[i][j] = (double **) malloc(sizeof(double) * ESTABLE_NRAD);

      for(int k=0; k<ESTABLE_NRAD; k++){
        dat->energy_data_storage[i][j][k] = (double *) malloc(sizeof(double) * ESTABLE_NENERGYBINS);
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
void free_extendedSourceTable(extendedSourceTable *tab) {
  if (tab != nullptr) {
    if (tab->dat != nullptr) {
      for (int ii = 0; ii < tab->n_a; ii++) {
        if (tab->dat[ii] != nullptr){
          free_extendedSourceData(tab->dat[ii]);
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
 * @details Dobby is free! (crossed out)
 * @param dat; extendedSourceData, data to be freed
 */
void free_extendedSourceData(extendedSourceData *dat) {
  if (dat != nullptr) {
    for (int ii =0; ii < ESTABLE_NHEIGHT; ii++) {
      if (dat->distance != nullptr) free(dat->distance[ii]);
      if (dat->flux != nullptr) free(dat->flux[ii]);
      if (dat->gshift != nullptr) free(dat->gshift[ii]);
      if (dat->refl_frac != nullptr) free(dat->refl_frac[ii]);
      if (dat->f_ad != nullptr) free(dat->f_ad[ii]);
      if (dat->f_inf != nullptr) free(dat->f_inf[ii]);
      if (dat->energy_data_storage != nullptr) free(dat->energy_data_storage[ii]);

      for (int jj = 0; jj < ESTABLE_NDISTANCE; jj++) {
        if (dat->flux[ii] != nullptr) free(dat->flux[ii][jj]);
        if (dat->gshift[ii] != nullptr) free(dat->gshift[ii][jj]);

        for (int kk = 0; kk < ESTABLE_NRAD; kk++){
          if (dat->energy_data_storage[ii][jj] != nullptr) free(dat->energy_data_storage[ii][jj][kk]);
        }
      }
    }
  }
    free(dat->radius);
    free(dat->height);
    free(dat->distance);
    free(dat->flux);
    free(dat->gshift);
    free(dat->energy_data_storage);
    free(dat->refl_frac);
    free(dat->f_inf);
    free(dat->f_ad);

}


/**
 *
 * @brief Calculates norm factor of the primary spectrum
 * @details Calculates normalisation factor of the primary spectrum to the number of photons reaching infinity
 * Need to add velocity of the source here
 * @param height; double, height
 * @param a; double, spin
 * @param gamma; double, photon index
 * @param f_inf; double, fraction of photons reaching infinity
 * @return double, normalisation factor
 */
double calc_norm_factor_primary_spectrum_ext(double height, double a, double gamma, double f_inf) {
  double g_inf = calc_g_inf(height, a);  // TODO: add velocity of the source (beta!)
  return f_inf / 0.5 * pow(g_inf, gamma);
}

/**
 *
 * Function which prints error status
 */
void printerror(int status) {
  printf("Error Code: %d\n",status);
}





