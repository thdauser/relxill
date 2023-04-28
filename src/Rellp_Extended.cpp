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

#include "cfitsio.h"
// #include "stdio.h"
#include "Relbase.h"
#include "Rellp_Extended.h"
#include "Rellp.h"

extern "C"{
#include "relutility.h"
#include "writeOutfiles.h"
}

extendedSourceTable *cached_extendedSource_table = NULL;

static extendedSourceTable* get_es_table(int* status){
  CHECK_STATUS_RET(*status,NULL);

  if (cached_extendedSource_table == NULL) {
    read_extendedSource_table(ESTABLE_FILENAME, &cached_extendedSource_table, status);
    CHECK_STATUS_RET(*status, NULL);
  }
  return cached_extendedSource_table;
}

static void get_extendedsourcetable_spin(int nrows, float **val, fitsfile *fptr, int *status){

  *val = (float*)  malloc(nrows* sizeof(float));

  for (int i=0; i<nrows; i++){

    int index = i + 2; // number of the row we read in the fits table
    fits_movabs_hdu(fptr, index, NULL, status);

    // dereference val to write a value to array
    fits_read_key(fptr, TFLOAT, "SPIN", &((*val)[i]), NULL, status);
    assert(*status == EXIT_SUCCESS);
  }

}

static extendedSourceData *
load_single_extendedSourceData(fitsfile *fptr, int n_h, int n_distance,
                               int n_rad, int n_storage, int *status,
                               double spin) {

  extendedSourceData *dat = new_extendedSourceData(status);
  CHECK_MALLOC_RET_STATUS(dat, status, NULL)


  int anynull = 0;
  double doublenull = 0.0;

  /*  read the columns */
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, n_rad, &doublenull, dat->radius,
                &anynull, status);

  for (int ii = 0; ii < n_h; ii++) {
    fits_read_col(fptr, TDOUBLE, 2, 1 + ii * n_h, 1, 1, &doublenull, &(dat->height[ii]),
                  &anynull, status);
    fits_read_col(fptr, TDOUBLE, 3, 1 + ii * n_h, 1, n_distance, &doublenull, dat->distance[ii],
                  &anynull, status);
  }

  for (int ii = 0; ii < n_h; ii++) {
    for (int jj = 0; jj < n_distance; jj++) {

      int irow = (ii * n_distance) + jj + 1;

      fits_read_col(fptr, TDOUBLE, 4, irow, 1, n_rad, &doublenull, dat->flux[ii][jj],
                    &anynull, status);

      for (int kk = 0; kk < n_rad; kk++){
        fits_read_col(fptr, TDOUBLE, 5, irow, kk*n_storage+1, n_storage, &doublenull, dat->energy_data_storage[ii][jj][kk],
                      &anynull, status);
      }

      // n_rad to 1 because it is not an array just a single value
      fits_read_col(fptr, TDOUBLE, 6, irow, 1, 1, &doublenull, &dat->f_ad[ii][jj],
                    &anynull, status);

      fits_read_col(fptr, TDOUBLE, 7, irow, 1, 1, &doublenull, &dat->f_inf[ii][jj],
                    &anynull, status);

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
void get_emis_extendedSource(emisProfile *emis, relParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  extendedSourceTable* estable = get_es_table(status);
  calc_extended_emis(emis, param, estable, status);

}

void read_extendedSource_table(char *filename, extendedSourceTable **inp_tab, int *status) {

  extendedSourceTable *tab = (*inp_tab);
  fitsfile *fptr = NULL;

  char fullfilename[999];

  do {
    if (tab != NULL) {
      RELXILL_ERROR("extended source table already loaded\n", status);
      break;
    }

    tab = new_extendedSourceTable(ESTABLE_NHEIGHT, ESTABLE_NRAD,
                                  ESTABLE_NDISTANCE, ESTABLE_NA, ESTABLE_NENERGYBINS, status);
    CHECK_STATUS_BREAK(*status);

    // should be set by previous routine
    assert(tab != NULL);
    assert(tab->dat != NULL);

    // get the full filename
    if (sprintf(fullfilename, "%s/%s", get_relxill_table_path(), filename) == -1) {
      RELXILL_ERROR("failed to construct full path the rel table\n", status);
      break;
    }

    // open the file
    if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
      CHECK_RELXILL_ERROR("opening of the extended source table failed", status);
      printf("    full path given: %s \n", filename);
      break;
    }

    // first read the spin axes of the table AND move to the correct extension
    printf("************************\n");
    printf("filename : %s\n", filename);


    get_extendedsourcetable_spin(tab->n_a, &(tab->a) , fptr, status);
    CHECK_RELXILL_ERROR("reading of spin axis failed\n", status);


    // now load the full table (need to go through all extensions
    int ii;
    int rownum;
    for (ii = 0; ii < tab->n_a; ii++) {

      rownum = ii + 2; // number of the row we read in the fits table
      fits_movabs_hdu(fptr, rownum, NULL, status);

      assert(tab->dat[ii] == NULL);

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

  if(fptr != NULL) {fits_close_file(fptr, status);}

}
void calc_extended_emis(emisProfile *emisProf, relParam *param,
                        extendedSourceTable *tab, int *status) {

  CHECK_STATUS_VOID(*status);

  emisProfile* emis_profile_table = interpol_estable(
      param->a, param->height, param->distance, param->gamma, ESTABLE_NENERGYBINS, tab, status);

  rebin_emisprofile_on_radial_grid(emisProf, emis_profile_table, status);

  free(emis_profile_table->re); // is not freed by free_emisProfile
  free_emisProfile(emis_profile_table);

  emisProf->returnFracs = calc_refl_frac_ext(param, tab, status);

  // take the beaming of the jet into account (see Dauser et al., 2013)
  // Maybe not necessary for extended source
//  for (int ii=0; ii<tab->n_rad; ii++){
//    if (param->beta > 1e-6) {
//      emisProf->emis[ii] *= pow(doppler_factor(emisProf->del_emit[ii], param->beta), 2);
//    }
//  }

  emisProf->normFactorPrimSpec =
      calc_norm_factor_primary_spectrum_ext(param->height, param->a, param->gamma, emisProf->returnFracs->f_inf);
}

emisProfile *interpol_estable(double a, double height, double distance,
                              double gamma, int n_gshift_storage,
                              extendedSourceTable *tab, int *status) {
  CHECK_STATUS_RET(*status, NULL);

  int ind_a;
  double ifac_a;
  get_ipol_factor((float) a, tab->a, tab->n_a, &ind_a, &ifac_a);

  extendedSourceData *dat_ind_a[2] = {tab->dat[ind_a], tab->dat[ind_a+1]};

  double* jet_rad = (double*) malloc(sizeof(double)*tab->n_rad);
  CHECK_MALLOC_RET_STATUS(jet_rad,status,NULL)
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
extendedSourceTable *new_extendedSourceTable(int n_h, int n_rad, int n_offaxis,
                                             int n_a, int n_storage,
                                             int *status) {
  extendedSourceTable *tab = (extendedSourceTable *) malloc(sizeof(extendedSourceTable));

  CHECK_MALLOC_RET_STATUS(tab, status, NULL)

  tab->n_a = n_a;
  tab->n_h = n_h;
  tab->n_rad = n_rad;
  tab->n_storage = n_storage;
  tab->n_offaxis = n_offaxis;

  tab->a = NULL;

  tab->dat = NULL;

  tab->dat = (extendedSourceData **) malloc(sizeof(extendedSourceData *) * tab->n_a);

  CHECK_MALLOC_RET_STATUS(tab->dat, status, tab)

  int ii;
  for (ii = 0; ii < tab->n_a; ii++) {
    tab->dat[ii] = NULL;
  }
  return tab;
}

extendedSourceData *new_extendedSourceData(int *status) {

  extendedSourceData *dat = (extendedSourceData *) malloc(sizeof(extendedSourceData));
  CHECK_MALLOC_RET_STATUS(dat, status, NULL)

  dat->height = (double *) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->height, status, NULL)

  dat->distance = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->distance, status, NULL)

  dat->radius = (double *) malloc(sizeof(double) * ESTABLE_NRAD);
  CHECK_MALLOC_RET_STATUS(dat->radius, status, NULL)

  dat->refl_frac = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, NULL)

  dat->f_ad = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->f_ad, status, NULL)

  dat->f_inf = (double **) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->f_inf, status, NULL)

  dat->flux = (double ***) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->flux, status, NULL)

  dat->gshift = (double ***) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->flux, status, NULL)

  dat->energy_data_storage = (double ****) malloc(sizeof(double) * ESTABLE_NHEIGHT);
  CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, NULL)

  for (int i=0; i<ESTABLE_NHEIGHT; i++){

    dat->refl_frac[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->refl_frac, status, NULL)

    dat->distance[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->distance, status, NULL)

    dat->f_ad[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->f_ad, status, NULL)

    dat->f_inf[i] = (double *) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->f_inf, status, NULL)

    dat->gshift[i] = (double **) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->gshift, status, NULL)

    dat->flux[i] = (double **) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->flux, status, NULL)

    dat->energy_data_storage[i] = (double ***) malloc(sizeof(double) * ESTABLE_NDISTANCE);
    CHECK_MALLOC_RET_STATUS(dat->energy_data_storage, status, NULL)

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

void free_extendedSourceTable(extendedSourceTable *tab) {
  if (tab != NULL) {
    if (tab->dat != NULL) {
      for (int ii = 0; ii < tab->n_a; ii++) {
        if (tab->dat[ii] != NULL){
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

void free_extendedSourceData(extendedSourceData *dat) {
  if (dat != NULL) {
    for (int ii =0; ii < ESTABLE_NHEIGHT; ii++) {
      if (dat->distance != NULL) free(dat->distance[ii]);
      if (dat->flux != NULL) free(dat->flux[ii]);
      if (dat->gshift != NULL) free(dat->gshift[ii]);
      if (dat->refl_frac != NULL) free(dat->refl_frac[ii]);
      if (dat->f_ad != NULL) free(dat->f_ad[ii]);
      if (dat->f_inf != NULL) free(dat->f_inf[ii]);
      if (dat->energy_data_storage != NULL) free(dat->energy_data_storage[ii]);

      for (int jj = 0; jj < ESTABLE_NDISTANCE; jj++) {
        if (dat->flux[ii] != NULL) free(dat->flux[ii][jj]);
        if (dat->gshift[ii] != NULL) free(dat->gshift[ii][jj]);

        for (int kk = 0; kk < ESTABLE_NRAD; kk++){
          if (dat->energy_data_storage[ii][jj] != NULL) free(dat->energy_data_storage[ii][jj][kk]);
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

double calc_norm_factor_primary_spectrum_ext(double height, double a, double gamma, double f_inf) {
  double g_inf = calc_g_inf(height, a);  // TODO: add velocity of the source (beta!)
  return f_inf / 0.5 * pow(g_inf, gamma);
}

void printerror(int status) {
  printf("Error Code: %d\n",status);
}





