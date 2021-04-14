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

#include "relreturn_table.h"

#include "xilltable.h"
#include "relbase.h"

returnTable *cached_retTable = NULL;

int global_rr_do_interpolation = 1;

/** create a new return table */
static returnTable *new_returnTable(int *status) {

  CHECK_STATUS_RET(*status, NULL);

  returnTable *tab = (returnTable *) malloc(sizeof(returnTable));
  CHECK_MALLOC_RET_STATUS(tab, status, tab)

  return tab;
}
/** init a new and empty return table (structure will be allocated)  */
static void init_returnTable(returnTable *tab, int nspin, int *status) {

  CHECK_STATUS_VOID(*status);

  tab->nspin = nspin;
  tab->spin = NULL;

  tab->retFrac = (tabulatedReturnFractions **) malloc(nspin * sizeof(tabulatedReturnFractions *));
  CHECK_MALLOC_VOID_STATUS(tab->retFrac, status)

}

void free_2d(double ***vals, int n1) {
  if (*vals != NULL) {
    for (int ii = 0; ii < n1; ii++) {
      free((*vals)[ii]);
    }
    free(*vals);
  }
}

static void free_returnFracData(tabulatedReturnFractions *dat) {

  if (dat != NULL) {

    free(dat->f_bh);
    free(dat->f_inf);
    free(dat->f_ret);

    free(dat->rlo);
    free(dat->rhi);

    free_2d(&(dat->frac_e), dat->nrad);
    free_2d(&(dat->frac_i), dat->nrad);

    free_2d(&(dat->gmin), dat->nrad);
    free_2d(&(dat->gmax), dat->nrad);

    if (dat->frac_g != NULL) {
      for (int ii = 0; ii < dat->nrad; ii++) {
        free_2d(&(dat->frac_g[ii]), dat->nrad);
      }
      free(dat->frac_g);
    }

    free(dat);
  }
}

static void free_returnTable(returnTable **tab) {

  if (*tab != NULL) {

    if ((*tab)->spin != NULL) {
      for (int ii = 0; ii < (*tab)->nspin; ii++) {
        free_returnFracData((*tab)->retFrac[ii]);
      }
      free((*tab)->spin);
      free((*tab)->retFrac);
    }
    free(*tab);
  }

  *tab = NULL;
}

void free_cached_returnTable(void) {
  free_returnTable(&cached_retTable);
}

static tabulatedReturnFractions *new_returnFracData(int nrad, int ng, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  tabulatedReturnFractions *dat = (tabulatedReturnFractions *) malloc(sizeof(tabulatedReturnFractions));
  CHECK_MALLOC_RET_STATUS(dat, status, dat)

  dat->nrad = nrad;
  dat->ng = ng;

  dat->f_bh = NULL;
  dat->f_ret = NULL;
  dat->f_inf = NULL;

  dat->frac_e = NULL;
  dat->frac_i = NULL;
  dat->frac_g = NULL;

  dat->gmin = NULL;
  dat->gmax = NULL;

  dat->rlo = NULL;
  dat->rhi = NULL;

  return dat;
}

void fits_moveToExtension(char *extname, fitsfile *fptr, int *status) {
  int extver = 0;
  fits_movnam_hdu(fptr, BINARY_TBL, extname, extver, status);
  if (*status != EXIT_SUCCESS) {
    printf(" *** error moving to extension %s\n", extname);
  }
}
/** read one axis of the rel table from the FITS file   */
static void get_reltable_axis_double(double **val,
                                     int *nval,
                                     char *extname,
                                     char *colname,
                                     fitsfile *fptr,
                                     int *status) {

  CHECK_STATUS_VOID(*status);

  fits_moveToExtension(extname, fptr, status);
  CHECK_STATUS_VOID(*status);


  // get the column id-number
  int colnum;
  if (fits_get_colnum(fptr, CASEINSEN, colname, &colnum, status)) return;

  long n;

  // get the number of rows
  if (fits_get_num_rows(fptr, &n, status)) return;

  // allocate memory for the array
  *val = (double *) malloc(n * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(*val, status)

  int anynul = 0;
  double nullval = 0.0;
  LONGLONG nelem = (LONGLONG) n;
  fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nelem, &nullval, *val, &anynul, status);

  *nval = (int) n;

}

static void get_returnRad_frac_dimensions(fitsfile *fptr, int *nr, int *ng, int *status) {

  CHECK_STATUS_VOID(*status);

  long n;

  // get the number of rows
  if (fits_get_num_rows(fptr, &n, status)) return;

  // read nr and ng from VARIABLES for now

  *nr = RETURNRAD_TABLE_NR;
  *ng = RETURNRAD_TABLE_NG;

  if (*nr != (int) n) {
    RELXILL_ERROR("return rad table: mismatch in number or radial bins ", status);
    printf("     expecting %i bins, but found %i rows in the table \n", *nr, (int) n);
  }

}

static void fits_rr_read_col(fitsfile *fptr,
                             double *val,
                             int firstrow,
                             int firstelem,
                             int nval,
                             int colnum,
                             int *status) {

  CHECK_STATUS_VOID(*status);

  int anynul = 0;
  double nullval = 0.0;
  LONGLONG nelem = (LONGLONG) nval;
  fits_read_col(fptr, TDOUBLE, colnum, (long) firstrow, (long) firstelem, nelem, &nullval, val, &anynul, status);

  relxill_check_fits_error(status);

}

int get_colnum(fitsfile *fptr, char *colname, int *status) {
  int colnum;
  fits_get_colnum(fptr, CASEINSEN, colname, &colnum, status);
  relxill_check_fits_error(status);
  return colnum;
}

static double *fits_rr_load_1d_data(fitsfile *fptr, char *colname, int nval, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  int colnum = get_colnum(fptr, colname, status);
  CHECK_STATUS_RET(*status, NULL);

  double *val = (double *) malloc(nval * sizeof(double));
  CHECK_MALLOC_RET_STATUS(val, status, NULL)

  fits_rr_read_col(fptr, val, 1, 1, nval, colnum, status);

  return val;
}

static double **fits_rr_load_2d_data(fitsfile *fptr, char *colname, int nval1, int nval2, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  int colnum = get_colnum(fptr, colname, status);
  CHECK_STATUS_RET(*status, NULL);

  double **val = (double **) malloc(nval1 * sizeof(double *));
  CHECK_MALLOC_RET_STATUS(val, status, NULL)

  for (int ii = 0; ii < nval1; ii++) {
    val[ii] = (double *) malloc(nval2 * sizeof(double));
    CHECK_MALLOC_RET_STATUS(val[ii], status, NULL)

    fits_rr_read_col(fptr, val[ii], ii + 1, 1, nval2, colnum, status);
  }

  return val;
}

static double ***fits_rr_load_3d_data(fitsfile *fptr, char *colname, int nval1, int nval2, int nval3, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  int colnum = get_colnum(fptr, colname, status);
  CHECK_STATUS_RET(*status, NULL);

  double ***val = (double ***) malloc(nval1 * sizeof(double **));
  CHECK_MALLOC_RET_STATUS(val, status, NULL)

  for (int ii = 0; ii < nval1; ii++) {
    val[ii] = (double **) malloc(nval2 * sizeof(double *));
    CHECK_MALLOC_RET_STATUS(val[ii], status, NULL)

    for (int jj = 0; jj < nval2; jj++) {
      val[ii][jj] = (double *) malloc(nval3 * sizeof(double));
      CHECK_MALLOC_RET_STATUS(val[ii][jj], status, NULL)

      fits_rr_read_col(fptr, val[ii][jj], ii + 1, jj * nval3 + 1, nval3, colnum, status);
    }
  }

  return val;
}

static tabulatedReturnFractions *fits_rr_load_single_fractions(fitsfile *fptr, char *extname, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  fits_moveToExtension(extname, fptr, status);
  CHECK_STATUS_RET(*status, NULL);

  int nrad = 0;
  int ng = 0;
  get_returnRad_frac_dimensions(fptr, &nrad, &ng, status);

  tabulatedReturnFractions *dat = new_returnFracData(nrad, ng, status);

  dat->rlo = fits_rr_load_1d_data(fptr, "rlo", nrad, status);
  dat->rhi = fits_rr_load_1d_data(fptr, "rhi", nrad, status);

  dat->frac_e = fits_rr_load_2d_data(fptr, "frac_e", nrad, nrad, status);
  dat->frac_i = fits_rr_load_2d_data(fptr, "frac_i", nrad, nrad, status);

  dat->gmin = fits_rr_load_2d_data(fptr, "gmin", nrad, nrad, status);
  dat->gmax = fits_rr_load_2d_data(fptr, "gmax", nrad, nrad, status);

  dat->frac_g = fits_rr_load_3d_data(fptr, "frac_g", nrad, nrad, ng, status);

  dat->f_ret = fits_rr_load_1d_data(fptr, "f_ret", nrad, status);
  dat->f_bh = fits_rr_load_1d_data(fptr, "f_ret", nrad, status);
  dat->f_inf = fits_rr_load_1d_data(fptr, "f_ret", nrad, status);

  return dat;
}

static void fits_rr_load_all_fractions(fitsfile *fptr, returnTable *tab, int *status) {

  CHECK_STATUS_VOID(*status);


  // currently our naming scheme only supports 99 spin values
  assert(tab->nspin <= 99);
  assert(tab->nspin > 0);
  char extname[50];

  for (int ii = 0; ii < tab->nspin; ii++) {
    sprintf(extname, "FRAC%02i", ii + 1);

    if (is_debug_run()) printf(" *** Return Rad: loading extension %s \n", extname);

    tab->retFrac[ii] = fits_rr_load_single_fractions(fptr, extname, status);
    CHECK_STATUS_BREAK(*status);

    tab->retFrac[ii]->a = tab->spin[ii]; // TODO: verify if we really need this
  }

}

static void fits_rr_load_returnRadTable(fitsfile *fptr, returnTable **inp_tab, int *status) {

  CHECK_STATUS_VOID(*status);

  returnTable *tab = new_returnTable(status);
  CHECK_STATUS_VOID(*status);

  /* get the number and values of the spin */
  double *spin;
  int nspin = 0;
  get_reltable_axis_double(&spin, &nspin, "SPIN", "a", fptr, status);

  /* initialize the table with them */
  init_returnTable(tab, nspin, status);
  tab->spin = spin;

  fits_rr_load_all_fractions(fptr, tab, status);

  (*inp_tab) = tab;

}

static void fits_read_returnRadTable(char *filename, returnTable **inp_tab, int *status) {

  CHECK_STATUS_VOID(*status);

  // open the table, stored at pwd or RELXILL_TABLE_PATH
  fitsfile *fptr = open_fits_table_stdpath(filename, status);

  // make sure we only store the table in a location which is empty / NULL
  assert(*inp_tab == NULL);

  fits_rr_load_returnRadTable(fptr, inp_tab, status);

  if (*status != EXIT_SUCCESS) {
    printf(" *** error *** initializing of the RETURN RADIATION table %s failed \n", filename);
    free_returnTable(inp_tab);
  }

  if (fptr != NULL) {
    fits_close_file(fptr, status);
  }

}

returnTable *get_returnrad_table(int *status) {

  if (cached_retTable==NULL) {
    fits_read_returnRadTable(RETURNRAD_TABLE_FILENAME, &cached_retTable, status);
  }

  return cached_retTable;
}


static int select_spinIndexForTable(double val_spin, double *arr_spin, int nspin, int *status) {

  // arr[k]<=val<arr[k+1]
  int k = binary_search(arr_spin, nspin, val_spin);

  // the spin of the table needs to be spin[k]>=val_spin
  if (arr_spin[k] < val_spin) {
    k++;
  }
  assert(arr_spin[k] >= val_spin);

  if (k >= nspin || k < 0) {
    RELXILL_ERROR("failed determining index of spin for the return radiation table \n", status);
    printf("    spin=%f leads to not allowed index of %i\n", val_spin, k);
  }

  return k;
}


static void allocate_radial_grid(returningFractions *ipol, double Rin, double Rout, int *status) {

  int klo_Rlo = binary_search(ipol->tabData->rlo, ipol->tabData->nrad, Rin);
  int khi_Rhi = binary_search(ipol->tabData->rhi, ipol->tabData->nrad, Rout);

  if( fabs(Rout - ipol->tabData->rhi[ipol->tabData->nrad-1]) < 1e-6){
    khi_Rhi=ipol->tabData->nrad-1;
  } else {
    khi_Rhi++;
  }

  if(fabs(Rin-ipol->tabData->rlo[0]) < 1e-6){
    klo_Rlo = 0;
  }


  int nrad_trim = (khi_Rhi+1) - klo_Rlo;

  assert(nrad_trim > 0);
  assert(nrad_trim <= ipol->tabData->nrad);

  ipol->nrad = nrad_trim;

  ipol->irad = (int *) malloc(nrad_trim * sizeof(int));
  CHECK_MALLOC_VOID_STATUS(ipol->irad, status)

  for (int ii = 0; ii < nrad_trim; ii++) {
    ipol->irad[ii] = klo_Rlo + ii;
  }

  ipol->rlo = (double *) malloc(nrad_trim * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(ipol->rlo, status)

  ipol->rhi = (double *) malloc(nrad_trim * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(ipol->rhi, status)

  ipol->rad = (double *) malloc(nrad_trim * sizeof(double));
  CHECK_MALLOC_VOID_STATUS(ipol->rad, status)

  ipol->proper_area_ring = NULL;

  for (int ii = 0; ii < nrad_trim; ii++) {
    ipol->rlo[ii] = ipol->tabData->rlo[ipol->irad[ii]];
    ipol->rhi[ii] = ipol->tabData->rhi[ipol->irad[ii]];
  }


  assert(Rin >= ipol->rlo[0]-1e-4); // TODO: make this limit stronger with a better table

  // reset lowest bin to Rin
  ipol->rlo[0] = Rin;
  ipol->rhi[nrad_trim-1] = Rout;

  for (int ii = 0; ii < nrad_trim; ii++) {
    ipol->rad[ii] = 0.5*(ipol->rlo[ii] + ipol->rhi[ii]);
  }


//  for (int ii = 0; ii < nrad_trim; ii++) {
//    ipol->proper_area_ring[ii] = calc_proper_area_ring(ipol->rlo[ii],ipol->rhi[ii], ipol->a);
//  }

}


static returningFractions *new_returningFractions(tabulatedReturnFractions *tab, double spin, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  returningFractions *dat = (returningFractions *) malloc(sizeof(returningFractions));
  CHECK_MALLOC_RET_STATUS(dat, status, dat)

  dat->tabData = tab;

  dat->a = spin;

  dat->irad = NULL;

  dat->proper_area_ring = NULL;

  dat->frac_i = NULL;
  dat->f_ret = NULL;
  dat->f_inf = NULL;

  return dat;
}

void free_returningFractions(returningFractions **dat) {

  if (*dat != NULL) {
    free((*dat)->rlo);
    free((*dat)->rhi);
    free((*dat)->rad);

    free_2d(&((*dat)->frac_i), (*dat)->nrad);

    free((*dat)->f_ret);
    free((*dat)->f_inf);

    free(*dat);
    *dat = NULL;
  }

}


static double **trim_fraci_to_radial_grid(double** tab_fraci, const int* ind_arr, int n, int *status) {

  double **fraci_trim = (double **) malloc(n * sizeof(double *));
  CHECK_MALLOC_RET_STATUS(fraci_trim, status, NULL)
  for (int ii = 0; ii < n; ii++) {
    fraci_trim[ii] = (double *) malloc(n * sizeof(double));
    CHECK_MALLOC_RET_STATUS(fraci_trim[ii], status, fraci_trim)

    for (int jj = 0; jj < n; jj++) {
      fraci_trim[ii][jj] = tab_fraci[ind_arr[ii]][ind_arr[jj]];
    }
  }


  return fraci_trim;
}


static double get_area_correction_factor(int index_radius, returningFractions* ret_fractions){
  double rlo_table = ret_fractions->tabData->rlo[ret_fractions->irad[index_radius]];

  if (ret_fractions->irad[index_radius]==0 && rlo_table>kerr_rms(ret_fractions->a)){ //TODO: will beremoved with dated table
    printf(" *** warning: resetting rlo to the ISCO at %f  (was %f before)\n",
           kerr_rms(ret_fractions->a), rlo_table);
    rlo_table = kerr_rms(ret_fractions->a);
  }


  double rhi_table = ret_fractions->tabData->rhi[ret_fractions->irad[index_radius]];
  double area_table = calc_proper_area_ring( rlo_table , rhi_table, ret_fractions->a);
  double area_model = calc_proper_area_ring( ret_fractions->rlo[index_radius] , ret_fractions->rhi[index_radius], ret_fractions->a);

  return area_model / area_table;
}


static double** get_interpolated_fraci(returningFractions *ret_fractions, int *status) {

  CHECK_STATUS_RET(*status,NULL);

  assert(ret_fractions->irad != NULL);
  double ** frac_i = trim_fraci_to_radial_grid(ret_fractions->tabData->frac_i,
                                               ret_fractions->irad,
                                               ret_fractions->nrad,
                                               status);

  int i_rad_rin = 0;
  int i_rad_rout = ret_fractions->nrad-1;

  double area_correction_rin = get_area_correction_factor(i_rad_rin, ret_fractions);
  double area_correction_rout = get_area_correction_factor(i_rad_rout, ret_fractions);

  assert( area_correction_rin-1  <1e-6);
  assert( area_correction_rout-1 <1e-6);


  if (global_rr_do_interpolation) {
    // correction of frac_e by emitted photons with respect to the smaller area of emission at i_rad_in
    for (int i_rad_incident = i_rad_rin; i_rad_incident < i_rad_rout; i_rad_incident++) {
      frac_i[i_rad_incident][i_rad_rin] *= area_correction_rin;
    }
    // same for i_rad_out
    for (int i_rad_incident = i_rad_rin; i_rad_incident < i_rad_rout-1; i_rad_incident++) {
      frac_i[i_rad_incident][i_rad_rout] *= area_correction_rout;
    }

  }

  return frac_i;
}


returningFractions *get_rrad_fractions(double spin, double rin, double rout, int *status) {

  CHECK_STATUS_RET(*status, NULL);
  assert(rin>0);
  assert(rout > rin);

  returnTable *tab = get_returnrad_table(status); // table will only be loaded if it is not already done

  int ind_spin = select_spinIndexForTable(spin, tab->spin, tab->nspin, status);
  tabulatedReturnFractions *tab_fractions = tab->retFrac[ind_spin];

  returningFractions* ret_fractions = new_returningFractions(tab_fractions, spin, status);


  allocate_radial_grid(ret_fractions, rin, rout, status);
  ret_fractions->frac_i = get_interpolated_fraci(ret_fractions, status);


  return ret_fractions;
}
