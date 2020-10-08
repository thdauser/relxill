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

#include "xilltable.h"

// possible parameters for the xillver tables
int global_param_index[] = {PARAM_GAM, PARAM_AFE, PARAM_LXI, PARAM_ECT, PARAM_KTE, PARAM_DNS,
                            PARAM_KTB, PARAM_ACO, PARAM_FRA, PARAM_INC};
char *global_param_names[] = {NAME_GAM, NAME_AFE, NAME_LXI, NAME_ECT, NAME_KTE, NAME_DNS,
                              NAME_KTB, NAME_ACO, NAME_FRA, NAME_INC};

// storage for the tables
xillTable *cached_xill_tab = NULL;
xillTable *cached_xill_tab_dens = NULL;
xillTable *cached_xill_tab_nthcomp = NULL;
xillTable *cached_xill_tab_dens_nthcomp = NULL;
xillTable *cached_xill_tab_ns = NULL;
xillTable *cached_xill_tab_co = NULL;

// todo: put "is model func" in relmodels.c and better structure?
static int is_ns_model(int model_type) {
  if ((model_type == MOD_TYPE_RELXILLNS) || (model_type == MOD_TYPE_XILLVERNS) || model_type == MOD_TYPE_RELXILLBBRET) {
    return 1;
  } else {
    return 0;
  }
}

static int is_co_model(int model_type) {
  if ((model_type == MOD_TYPE_RELXILLCO) || (model_type == MOD_TYPE_XILLVERCO)) {
    return 1;
  } else {
    return 0;
  }
}

static int is_CpD_model(int model_type) {
  if ((model_type == MOD_TYPE_RELXILLDENS_NTHCOMP) ||
      (model_type == MOD_TYPE_RELXILLLPDENS_NTHCOMP) ||
      (model_type == MOD_TYPE_XILLVERDENS_NTHCOMP)
      ) {
    return 1;
  } else {
    return 0;
  }
}

static int get_num_elem(const int *n_parvals, int npar) {

  assert(npar >= 1);
  long num_elem = n_parvals[0];

  int ii;
  for (ii = 1; ii < npar; ii++) {
    num_elem *= n_parvals[ii];
  }

  assert(num_elem <= INT_MAX);

  return (int) num_elem;
}

static void init_xilltable_data_struct(xillTable *tab, int *status) {

  CHECK_STATUS_VOID(*status);

  tab->num_elements = get_num_elem(tab->num_param_vals, tab->num_param);

  tab->data_storage = (float **) malloc(sizeof(float *) * tab->num_elements);
  CHECK_MALLOC_VOID_STATUS(tab->data_storage, status)

  // important to make sure everything is set to NULL (used to only load spectra if !=NULL)
  int ii;
  for (ii = 0; ii < tab->num_elements; ii++) {
    tab->data_storage[ii] = NULL;
  }

}

/** get a new and empty rel table (structure will be allocated)  */
xillTable *new_xillTable(int num_param, int *status) {

  xillTable *tab = (xillTable *) malloc(sizeof(xillTable));
  CHECK_MALLOC_RET_STATUS(tab, status, tab)

  tab->num_param = num_param;

  // first make sure we ask only for the dimensions that are implemented
  if (!(tab->num_param == 5 || tab->num_param == 6)) {
    RELXILL_ERROR("wrong dimensionality of the xillver table", status);
    return NULL;
  }

  tab->num_param_vals = (int *) malloc(sizeof(int) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->num_param_vals, status, NULL)

  tab->param_index = (int *) malloc(sizeof(int) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->param_index, status, NULL)

  tab->param_names = (char **) malloc(sizeof(char *) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->param_index, status, NULL)

  tab->param_vals = (float **) malloc(sizeof(float *) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->param_vals, status, NULL)

  int ii;
  for (ii = 0; ii < num_param; ii++) {
    tab->param_vals[ii] = NULL;
    tab->param_names[ii] = NULL;
  }

  tab->elo = NULL;
  tab->ehi = NULL;

  // inclination is the last value
  tab->incl = NULL;
  tab->n_incl = -1;

  tab->data_storage = NULL;

  return tab;
}

static int is_dens_model(int model_type) {
  if ((model_type == MOD_TYPE_RELXILLDENS) || (model_type == MOD_TYPE_XILLVERDENS) ||
      (model_type == MOD_TYPE_RELXILLLPDENS)) {
    return 1;
  } else {
    return 0;
  }
}

int is_6dim_table(int model_type) {
  if (is_co_model(model_type)) {
    return 1;
  } else {
    return 0;
  }
}

static int get_num_param_auto(fitsfile *fptr, int *status) {

  CHECK_STATUS_RET(*status, 0);

  // get the number of rows
  long n;
  fits_get_num_rows(fptr, &n, status);
  relxill_check_fits_error(status);

  return (int) n;
}

/* routine to get the pointer to the pre-defined parameter index array
 * (i.e. with this we know which parameter relates to the input parameters ) */
static void set_parindex_from_parname(int *pindex, char **pname, int n, int *status) {

  int ii;
  // loop over given parameters
  for (ii = 0; ii < n; ii++) {
    pindex[ii] = -1;

    // now loop over all possible parameters
    int jj;
    for (jj = 0; jj < N_PARAM_MAX; jj++) {
      if (strcmp(pname[ii], global_param_names[jj]) == 0) {
        pindex[ii] = global_param_index[jj];
        continue;
      }
    }

    if (pindex[ii] == -1) {
      RELXILL_ERROR(" parameter given in xillver table not found\n", status);
      printf("    parameter ** %s ** from xillver table, not known to relxill \n", pname[ii]);
      printf("    please make sure you downloaded the correct table \n");
    }
  }

}

void print_xilltable_parameters(const xillTable *tab, char *const *xilltab_parname) {
  int ii;
  for (ii = 0; ii < tab->num_param; ii++) {
    printf(" loaded parameter %s  (index=%i) \t -  %02i values from %.2f to %.2f\n",
           xilltab_parname[ii], tab->param_index[ii], tab->num_param_vals[ii],
           tab->param_vals[ii][0], tab->param_vals[ii][tab->num_param_vals[ii] - 1]);
  }
}

/** read the parameters of the xillver FITS table   */
static void get_xilltable_parameters(fitsfile *fptr, xillTable *tab, xillParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  int extver = 0;
  fits_movnam_hdu(fptr, BINARY_TBL, "PARAMETERS", extver, status);
  if (*status != EXIT_SUCCESS) {
    printf(" *** error moving to extension PARAMETERS in the xillver table\n");
    return;
  }

  // we know the column for the Number of Parameter-Values
  int colnum_n = 9;
  int colnum_vals = 10;

  // get the number of rows
  long n;
  if (fits_get_num_rows(fptr, &n, status)) return;

  if (tab->num_param != n) {
    RELXILL_ERROR("wrong format of the xillver table (not enough or too many parameter values tabulated)", status);
  }

  int anynul = 0;
  double nullval = 0.0;

  int ii;
  char strnull[10];
  strcpy(strnull, " ");

  /* allocate space for string column value */
  char **xilltab_parname = (char **) malloc(sizeof(char *) * tab->num_param);
  CHECK_MALLOC_VOID_STATUS(xilltab_parname, status)

  // get the name
  for (ii = 0; ii < tab->num_param; ii++) {
    tab->param_names[ii] = (char *) malloc(sizeof(char) * 8);
    CHECK_MALLOC_VOID_STATUS(tab->param_names, status)
  }
  fits_read_col(fptr, TINT, colnum_n, 1, 1, tab->num_param, &nullval, tab->num_param_vals, &anynul, status);
  fits_read_col(fptr, TSTRING, 1, 1, 1, tab->num_param, strnull, tab->param_names, &anynul, status);

  for (ii = 0; ii < tab->num_param; ii++) {
    /** the we load the parameter values **/
    float **ptr_val = &(tab->param_vals[ii]);

    if (ptr_val == NULL) {
      RELXILL_ERROR(" *** error loading the xillver parameter from the table", status);
      return;
    }
    (*ptr_val) = (float *) malloc(tab->num_param_vals[ii] * sizeof(float));
    CHECK_MALLOC_VOID_STATUS(*ptr_val, status)

    fits_read_col(fptr, TFLOAT, colnum_vals, ii + 1, 1, tab->num_param_vals[ii], &nullval, *ptr_val, &anynul,
                  status);

    CHECK_STATUS_BREAK(*status);
  }

  // set the index of each parameter according to the NAME column we found in the table
  set_parindex_from_parname(tab->param_index, tab->param_names, tab->num_param, status);

  if (is_debug_run()) {
    print_xilltable_parameters(tab, tab->param_names);
  }

  // now we set a pointer to the inclination separately (assuming it is the last parameter)
  tab->incl = tab->param_vals[tab->num_param - 1];
  tab->n_incl = tab->num_param_vals[tab->num_param - 1];

}

// static void get_reltable_axis(int nrows, float** val, char* extname, char* colname, fitsfile* fptr, int* status){
static void get_xilltable_ener(int *n_ener, float **elo, float **ehi, fitsfile *fptr, int *status) {

  CHECK_STATUS_VOID(*status);

  int extver = 0;
  fits_movnam_hdu(fptr, BINARY_TBL, "ENERGIES", extver, status);
  if (*status != EXIT_SUCCESS) {
    printf(" *** error moving to extension ENERGIES\n");
    return;
  }

  // get the column id-number
  int colnum_elo = 1;
  int colnum_ehi = 2;

  // get the number of rows
  long nrows;
  if (fits_get_num_rows(fptr, &nrows, status)) return;
  // (strongly assume they fit in Integer range)
  *n_ener = (int) nrows;

  // allocate memory for the array
  *elo = (float *) malloc((*n_ener) * sizeof(float));
  CHECK_MALLOC_VOID_STATUS(*elo, status)
  *ehi = (float *) malloc((*n_ener) * sizeof(float));
  CHECK_MALLOC_VOID_STATUS(*ehi, status)

  int anynul = 0;
  double nullval = 0.0;
  LONGLONG nelem = (LONGLONG) (*n_ener);
  fits_read_col(fptr, TFLOAT, colnum_elo, 1, 1, nelem, &nullval, *elo, &anynul, status);
  fits_read_col(fptr, TFLOAT, colnum_ehi, 1, 1, nelem, &nullval, *ehi, &anynul, status);

  relxill_check_fits_error(status);

  // simply let's check if something went wrong
  CHECK_RELXILL_ERROR("reading of energy grid of the xillver table failed", status);

}

static float *get_xill_param_vals_array(xillParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  float *param_vals = (float *) malloc(N_PARAM_MAX * sizeof(float));
  CHECK_MALLOC_RET_STATUS(param_vals, status, NULL)
  // store the input in a variable
  param_vals[PARAM_GAM] = (float) param->gam;
  param_vals[PARAM_AFE] = (float) param->afe;
  param_vals[PARAM_LXI] = (float) param->lxi;
  param_vals[PARAM_ECT] = (float) param->ect;
  param_vals[PARAM_DNS] = (float) param->dens;
  param_vals[PARAM_KTB] = (float) param->kTbb;
  param_vals[PARAM_FRA] = (float) param->frac_pl_bb;
  param_vals[PARAM_INC] = (float) param->incl;

  return param_vals;
}

static int *xillInd_from_parInput(xillParam *param, xillTable *tab, int *status) {

  CHECK_STATUS_RET(*status, NULL);
  // store the input in a variable
  float *inp_param_vals = get_xill_param_vals_array(param, status);

  int ii;
  int *ind = (int *) malloc(XILLTABLE_N_PARAM * sizeof(int));
  CHECK_MALLOC_RET_STATUS(ind, status, NULL)

  for (ii = 0; ii < tab->num_param; ii++) {

    ind[ii] = binary_search_float(tab->param_vals[ii], tab->num_param_vals[ii],
                                  inp_param_vals[tab->param_index[ii]]);

    // make sure all parameters are by default within the defined limits here!!
    if (ind[ii] < 0) {
      ind[ii] = 0;
    } else if (ind[ii] > tab->num_param_vals[ii] - 2) {
      ind[ii] = tab->num_param_vals[ii] - 2;
    }
  }

  free(inp_param_vals);

  return ind;
}

static int
get_xillspec_rownum(const int *num_param_vals, int num_param, int nn, int ii, int jj, int kk, int ll, int mm) {

  if (num_param == 5) {
    return (((ii * num_param_vals[1] + jj) * num_param_vals[2]
        + kk) * num_param_vals[3] + ll) * num_param_vals[4] + mm + 1;
  } else if (num_param == 6) {
    return ((((nn * num_param_vals[1] + ii) * num_param_vals[2]
        + jj) * num_param_vals[3] + kk) * num_param_vals[4] + ll) * num_param_vals[5] + mm + 1;
  }

  return -1;
}

static void set_dat(float *spec, xillTable *tab, int i0, int i1, int i2, int i3, int i4, int i5) {
  int index = get_xillspec_rownum(tab->num_param_vals, tab->num_param,
                                  i0, i1, i2, i3, i4, i5);
  tab->data_storage[index] = spec;
}

// get one Spectrum from the Data Storage
static float *get_dat(xillTable *tab, int i0, int i1, int i2, int i3, int i4, int i5) {

  int index = get_xillspec_rownum(tab->num_param_vals, tab->num_param,
                                  i0, i1, i2, i3, i4, i5);

  return tab->data_storage[index];
}

char *getFullPathTableName(char *filename, int *status) {

  int MAXSIZE = 1000;
  char *fullfilename = malloc(sizeof(char *) * MAXSIZE);

  if (sprintf(fullfilename, "%s/%s", get_relxill_table_path(), filename) == -1) {
    RELXILL_ERROR("failed to construct full path the rel table", status);
    return NULL;
  }

  return fullfilename;

}

int checkIfTableExists(char *filename, int *status) {

  char *fullfilename = getFullPathTableName(filename, status);
  CHECK_STATUS_RET(*status, 0);

  fitsfile *fptr = NULL;
  int statusTableExists = EXIT_SUCCESS;
  int tableExists = 0;

  if (fits_open_table(&fptr, fullfilename, READONLY, &statusTableExists) == 0) {
    tableExists = 1;
  }

  if (fptr != NULL) {
    fits_close_file(fptr, &statusTableExists);
  }

  return tableExists;
}

fitsfile *open_fits_table_stdpath(char *filename, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  char *fullfilename = getFullPathTableName(filename, status);
  CHECK_STATUS_RET(*status, NULL);

  fitsfile *fptr = NULL;
  if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
    RELXILL_ERROR("opening of the table failed", status);
    printf("    either the full path given (%s) is wrong \n", fullfilename);
    printf("    or you need to download the table ** %s **  from \n", filename);
    printf("    http://www.sternwarte.uni-erlangen.de/research/relxill/ \n");
  }

  free(fullfilename);

  return fptr;
}

/** load the complete relline table */
void init_xillver_table(char *filename, xillTable **inp_tab, xillParam *param, int *status) {

  CHECK_STATUS_VOID(*status);

  xillTable *tab = (*inp_tab);
  fitsfile *fptr = NULL;

  print_version_number(status);

  fptr = open_fits_table_stdpath(filename, status);
  CHECK_STATUS_VOID(*status);

  assert(tab == NULL);

  int num_param = get_num_param_auto(fptr, status);
  CHECK_STATUS_VOID(*status);

  /** =1= allocate space for the new table  **/
  tab = new_xillTable(num_param, status);

  /** =2= now load the energy grid **/
  get_xilltable_ener(&(tab->n_ener), &(tab->elo), &(tab->ehi), fptr, status);

  /** =3= and now the stored parameter values (also check if the correct number of parameters) **/
  get_xilltable_parameters(fptr, tab, param, status);

  /** =4= finally set up the complete data structure **/
  init_xilltable_data_struct(tab, status);

  if (*status == EXIT_SUCCESS) {
    // assign the value
    (*inp_tab) = tab;
  } else {
    printf(" *** error *** initializing of the XILLVER table %s failed \n", filename);
    free_xillTable(tab);
  }

  if (fptr != NULL) {
    fits_close_file(fptr, status);
  }

}

/** check if the given parameter index is in the xillver table
 *  -1 : not found
 *  ii : return index where this value is found in the param_index array  **/
static int getParameterTableIndex(xillTable *tab, int ind) {

  int ii;
  for (ii = 0; ii < tab->num_param; ii++) {
    if (tab->param_index[ii] == ind) {
      return ii;
    }
  }
  return -1;
}

// convert the index from the (actual) xillver FITS table parameter
// extension to out (strange!) universal parameter index which can
// handle 5dim and 6dim tables
static int convertTo6dimTableIndex(int num_param, int indexTable) {
  assert(num_param == 6 || num_param == 5);
  if (num_param == 5) {
    indexTable++; // for 5dim, the first index is left empty
  }
  return indexTable;
}

static void normalizeXillverSpecLogxiDensity(float *spec,
                                             xillTable *tab,
                                             double density,
                                             double logxi,
                                             const int *indVals,
                                             int array) {

  int ind_dens = getParameterTableIndex(tab, PARAM_DNS);
  int ind_lxi = getParameterTableIndex(tab, PARAM_LXI);

  int ind_dens_6dim = convertTo6dimTableIndex(tab->num_param, ind_dens);
  int ind_lxi_6dim = convertTo6dimTableIndex(tab->num_param, ind_lxi);

  // ind==-1 if parameter not contained
  if (ind_dens >= 0) {
    density = tab->param_vals[ind_dens][indVals[ind_dens_6dim]];

  }
  if (ind_lxi >= 0) {
    logxi = tab->param_vals[ind_lxi][indVals[ind_lxi_6dim]];
  }

  renorm_xill_spec(spec, tab->n_ener, logxi, density);

}

static void load_single_spec(char *fname, fitsfile **fptr, xillTable *tab, double defDensity, double defLogxi,
                             int nn, int ii, int jj, int kk, int ll, int mm, int *status) {

  CHECK_STATUS_VOID(*status);

  // open the fits file if not already open
  if (*fptr == NULL) {
    *fptr = open_fits_table_stdpath(fname, status);
    CHECK_STATUS_VOID(*status);
  }

  int extver = 0;
  fits_movnam_hdu(*fptr, BINARY_TBL, "SPECTRA", extver, status);
  if (*status != EXIT_SUCCESS) {
    printf(" *** error moving to extension SPECTRA in the xillver table\n");
    return;
  }

  float *spec = (float *) malloc(tab->n_ener * sizeof(float));
  CHECK_MALLOC_VOID_STATUS(spec, status)

  int colnum_spec = 2;
  int anynul = 0;
  double nullval = 0.0;
  LONGLONG nelem = (LONGLONG) tab->n_ener;

  // get the row number (this is how the xillver table is defined)
  int rownum = get_xillspec_rownum(tab->num_param_vals, tab->num_param, nn, ii, jj, kk, ll, mm);

  fits_read_col(*fptr, TFLOAT, colnum_spec, rownum, 1, nelem, &nullval, spec, &anynul, status);
  if (*status != EXIT_SUCCESS) {
    printf("\n *** ERROR *** failed reading table %s  (rownum %i) \n",
           fname, rownum);
    relxill_check_fits_error(status);
    return;
  }

  const int n_array = 6;
  int indexArray[] = {nn, ii, jj, kk, ll, mm};
  normalizeXillverSpecLogxiDensity(spec, tab, defDensity, defLogxi, indexArray, n_array);

  set_dat(spec, tab, nn, ii, jj, kk, ll, mm);
}

// Default values from the param-structure, as those are set to the value the table is calculated
// with, even if the model does not have a lxi or density parameter
static double getDefaultDensity(xillParam *param) {
  return param->dens;
}
static double getDefaultLogxi(xillParam *param) {
  return param->lxi;
}

static void check_xillTable_cache(char *fname, xillParam *param, xillTable *tab, const int *ind, int *status) {

  CHECK_STATUS_VOID(*status);

  // =2=  check if the necessary spectra are loaded (we only open the file once)
  fitsfile *fptr = NULL;
  int ii;
  int jj;
  int kk;
  int ll;
  int mm;
  int nn;

  // (current) standard case for 5 param
  int i0lo = 0;
  int i0hi = 0;
  int istart = 0;
  // for 6dim we (i.e., I) decided to use the complicated statement
  // (reason: new dimension has to be at the beginning, as inclination HAS to be last
  if (tab->num_param == 6) {
    i0lo = ind[0];
    i0hi = ind[0] + 1;
    istart = 1;
  }

  double defDensity = getDefaultDensity(param);
  double defLogxi = getDefaultLogxi(param);

  for (nn = i0lo; nn <= i0hi; nn++) { // for 5dim this is a dummy loop
    for (ii = ind[istart]; ii <= ind[istart] + 1; ii++) {
      for (jj = ind[istart + 1]; jj <= ind[istart + 1] + 1; jj++) {
        for (kk = ind[istart + 2]; kk <= ind[istart + 2] + 1; kk++) {
          for (ll = ind[istart + 3]; ll <= ind[istart + 3] + 1; ll++) {
            // always load **all** incl bins as for relxill we will certainly need it
            for (mm = 0; mm < tab->n_incl; mm++) {

              if (get_dat(tab, nn, ii, jj, kk, ll, mm) == NULL) {
                load_single_spec(fname, &fptr, tab, defDensity, defLogxi, nn, ii, jj, kk, ll, mm, status);
                CHECK_STATUS_VOID(*status);
              }

            }
          }
        }
      }
    }

  }

  if (fptr != NULL) {
    if (fits_close_file(fptr, status)) {
      RELXILL_ERROR(" *** error closing FITS file", status);
    }
  }

}

xill_spec *new_xill_spec(int n_incl, int n_ener, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  xill_spec *spec = (xill_spec *) malloc(sizeof(xill_spec));
  CHECK_MALLOC_RET_STATUS(spec, status, NULL)

  spec->n_ener = n_ener;
  spec->n_incl = n_incl;

  spec->ener = (double *) malloc(sizeof(double) * (n_ener + 1));
  CHECK_MALLOC_RET_STATUS(spec->ener, status, spec)

  spec->incl = (double *) malloc(sizeof(double) * (n_incl));
  CHECK_MALLOC_RET_STATUS(spec->incl, status, spec)

  spec->flu = (double **) malloc(sizeof(double *) * n_incl);
  CHECK_MALLOC_RET_STATUS(spec->flu, status, spec)

  int ii;
  for (ii = 0; ii < n_incl; ii++) {
    spec->flu[ii] = (double *) malloc(sizeof(double) * n_ener);
    CHECK_MALLOC_RET_STATUS(spec->flu[ii], status, spec)
  }

  assert(spec != NULL);

  return spec;
}

void free_xill_spec(xill_spec *spec) {
  if (spec != NULL) {
    free(spec->ener);
    free(spec->incl);
    if (spec->flu != NULL) {
      int ii;
      for (ii = 0; ii < spec->n_incl; ii++) {
        free(spec->flu[ii]);
      }
      free(spec->flu);
    }
    free(spec);
  }
}

static void interp_5d_tab_incl(xillTable *tab, double *flu, int n_ener,
                               double f1, double f2, double f3, double f4,
                               int i0, int i1, int i2, int i3, int i4, int i5) {

  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    flu[ii] =
        ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2, i3, i4, i5))[ii]) +
            (f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3, i4, i5))[ii]) +
            (1.0 - f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3, i4, i5))[ii]) +
            (1.0 - f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2, i3 + 1, i4, i5))[ii]) +
            (f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3, i4, i5))[ii]) +
            (f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3 + 1, i4, i5))[ii]) +
            (1.0 - f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3 + 1, i4, i5))[ii]) +
            (f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4, i5))[ii]))
            * (1 - f4) +
            ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2, i3, i4 + 1, i5))[ii]) +
                (f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3, i4 + 1, i5))[ii]) +
                (1.0 - f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3, i4 + 1, i5))[ii]) +
                (1.0 - f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2, i3 + 1, i4 + 1, i5))[ii]) +
                (f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3, i4 + 1, i5))[ii]) +
                (f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3 + 1, i4 + 1, i5))[ii]) +
                (1.0 - f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3 + 1, i4 + 1, i5))[ii]) +
                (f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4 + 1, i5))[ii]))
                * (f4);
  }

}

static void interp_5d_tab(xillTable *tab, double *flu, int n_ener,
                          double f1, double f2, double f3, double f4, double f5,
                          int i0, int i1, int i2, int i3, int i4,
                          int i5) {

  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    flu[ii] =
        (((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2, i3, i4, i5))[ii]) +
            (f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3, i4, i5))[ii]) +
            (1.0 - f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3, i4, i5))[ii]) +
            (1.0 - f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2, i3 + 1, i4, i5))[ii]) +
            (f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3, i4, i5))[ii]) +
            (f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3 + 1, i4, i5))[ii]) +
            (1.0 - f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3 + 1, i4, i5))[ii]) +
            (f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4, i5))[ii]))
            * (1 - f4) +
            ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2, i3, i4 + 1, i5))[ii]) +
                (f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3, i4 + 1, i5))[ii]) +
                (1.0 - f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3, i4 + 1, i5))[ii]) +
                (1.0 - f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2, i3 + 1, i4 + 1, i5))[ii]) +
                (f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3, i4 + 1, i5))[ii]) +
                (f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3 + 1, i4 + 1, i5))[ii]) +
                (1.0 - f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3 + 1, i4 + 1, i5))[ii]) +
                (f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4 + 1, i5))[ii]))
                * (f4)) * (1.0 - f5) +
            (((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2, i3, i4, i5 + 1))[ii]) +
                (f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3, i4, i5 + 1))[ii]) +
                (1.0 - f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3, i4, i5 + 1))[ii]) +
                (1.0 - f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2, i3 + 1, i4, i5 + 1))[ii]) +
                (f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3, i4, i5 + 1))[ii]) +
                (f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3 + 1, i4, i5 + 1))[ii]) +
                (1.0 - f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3 + 1, i4, i5 + 1))[ii]) +
                (f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4, i5 + 1))[ii]))
                * (1 - f4) +
                ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2, i3, i4 + 1, i5 + 1))[ii]) +
                    (f1) * (1.0 - f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3, i4 + 1, i5 + 1))[ii]) +
                    (1.0 - f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3, i4 + 1, i5 + 1))[ii]) +
                    (1.0 - f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2, i3 + 1, i4 + 1, i5 + 1))[ii]) +
                    (f1) * (f2) * (1.0 - f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3, i4 + 1, i5 + 1))[ii]) +
                    (f1) * (1.0 - f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2, i3 + 1, i4 + 1, i5 + 1))[ii]) +
                    (1.0 - f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1, i2 + 1, i3 + 1, i4 + 1, i5 + 1))[ii]) +
                    (f1) * (f2) * (f3) * (double) ((get_dat(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4 + 1, i5 + 1))[ii]))
                    * (f4)) * (f5);

  }

}

/* cheap man's version of a 6dim interpolation */
static void interp_6d_tab_incl(xillTable *tab, double *flu, int n_ener,
                               double *fac, int nfac, const int *ind, int iincl) {

  assert(nfac == 6);

  double s1[n_ener];
  double s2[n_ener];

  interp_5d_tab_incl(tab, s1, n_ener,
                     fac[1], fac[2], fac[3], fac[4], ind[0],
                     ind[1], ind[2], ind[3], ind[4], iincl);

  interp_5d_tab_incl(tab, s2, n_ener,
                     fac[1], fac[2], fac[3], fac[4], ind[0] + 1,
                     ind[1], ind[2], ind[3], ind[4], iincl);

  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    flu[ii] = interp_lin_1d(fac[0], s1[ii], s2[ii]);
  }
}

/* cheap man's version of a 6dim interpolation */
static void interp_6d_tab(xillTable *tab, double *flu, int n_ener,
                          double *fac, int nfac, const int *ind) {

  assert(nfac == 6);

  double s1[n_ener];
  double s2[n_ener];

  interp_5d_tab(tab, s1, n_ener,
                fac[1], fac[2], fac[3], fac[4], fac[5],
                ind[0],
                ind[1], ind[2], ind[3], ind[4], ind[5]);

  interp_5d_tab(tab, s2, n_ener,
                fac[1], fac[2], fac[3], fac[4], fac[5],
                ind[0] + 1,
                ind[1], ind[2], ind[3], ind[4], ind[5]);

  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    flu[ii] = interp_lin_1d(fac[0], s1[ii], s2[ii]);
  }
}

/** check boundary of the Ecut parameter  (grav. redshift can lead to the code asking
 *  for an Ecut not tabulated, although actually observed ecut is larger
 *  Input: tab, param
 *  Output fac         **/
static void check_boundarys_ecut(xillTable *tab, const xillParam *param, double *fac) {
  // first need to make sure ECUT is a parameter in the table
  int index_ect = getParameterTableIndex(tab, PARAM_ECT);
  if (index_ect >= 0) {
    if (param->ect <= tab->param_vals[tab->param_index[index_ect]][0]) {
      fac[index_ect] = 0.0;
    }
    // can happen due to grav. redshift, although actually observed ecut is larger
    if (param->ect >=
        tab->param_vals[tab->param_index[index_ect]][tab->num_param_vals[tab->param_index[index_ect]] - 1]) {
      fac[index_ect] = 1.0;
    }
  }
}

static void resetInpvalsToBoundaries(char *pname, float *p_inpVal, float tabValLo, float tabValHi) {

  float inpVal = *p_inpVal;
  if (inpVal < tabValLo) {
    if (is_debug_run()) {
      printf(" *** warning: paramter %s=%e below lowest table value, restetting to %e\n",
             pname, inpVal, tabValLo);
    }
    inpVal = tabValLo;
  } else if (inpVal > tabValHi) {
    if (is_debug_run()) {
      printf("\n *** warning: paramter %s=%e above largest table value, restetting to %e\n",
             pname, inpVal, tabValHi);
    }
    inpVal = tabValHi;
  }

  *p_inpVal = inpVal;

}

static xill_spec *interp_xill_table(xillTable *tab, xillParam *param, const int *ind, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  xill_spec *spec = NULL;
  if (is_xill_model(param->model_type)) {
    spec = new_xill_spec(1, tab->n_ener, status);
  } else {
    spec = new_xill_spec(tab->n_incl, tab->n_ener, status);
  }

  assert(spec != NULL);
  assert(spec->n_ener == tab->n_ener);

  // set the energy grid
  int ii;
  for (ii = 0; ii < spec->n_ener; ii++) {
    spec->ener[ii] = tab->elo[ii];
  }
  spec->ener[spec->n_ener] = tab->ehi[spec->n_ener - 1];

  // set the inclination grid
  for (ii = 0; ii < spec->n_incl; ii++) {
    spec->incl[ii] = tab->incl[ii];
  }

  float *inp_param_vals = get_xill_param_vals_array(param, status);

  int nfac = tab->num_param;
  double *fac = (double *) malloc(sizeof(double) * nfac);
  CHECK_MALLOC_RET_STATUS(fac, status, NULL)

  /* calculate interpolation factor for all parameters
   * ([nfac-1] is inclination, which might not be used ) */
  int pind;
  for (ii = 0; ii < nfac; ii++) {
    // need the index
    pind = tab->param_index[ii];

    if ( ! (( !is_xill_model(param->model_type)) &&  (tab->param_index[ii] == PARAM_INC ))) {
      resetInpvalsToBoundaries(tab->param_names[ii], &(inp_param_vals[pind]), tab->param_vals[ii][0],
                               tab->param_vals[ii][tab->num_param_vals[ii] - 1]);
    }

    fac[ii] = (inp_param_vals[pind] - tab->param_vals[ii][ind[ii]]) /
        (tab->param_vals[ii][ind[ii] + 1] - tab->param_vals[ii][ind[ii]]);

    if (is_debug_run()) {
      printf("\n [%i] %s (par_index=%i)  : parval = %.2e (index=%i), determining fac = %.2e ",
             ii, tab->param_names[ii], pind, inp_param_vals[pind], ind[ii], fac[ii]);
    }
  }

  free(inp_param_vals);

  /** check boundary of the Ecut parameter
   *  (grav. redshift can lead to the code asking for an Ecut not tabulated,
   *  although actually observed ecut is larger **/
  check_boundarys_ecut(tab, param, fac);

  if (tab->num_param == 5) {
    if (is_xill_model(param->model_type)) {
      interp_5d_tab(tab, spec->flu[0], spec->n_ener,
                    fac[0], fac[1], fac[2], fac[3], fac[4], 0,
                    ind[0], ind[1], ind[2], ind[3], ind[4]);

    } else {
      // do not interpolate over the inclination (so skip the last parameter)
      nfac--;
      assert(nfac == 4);
      // get the spectrum for EACH flux bin
      for (ii = 0; ii < spec->n_incl; ii++) {
        interp_5d_tab_incl(tab, spec->flu[ii], spec->n_ener,
                           fac[0], fac[1], fac[2], fac[3], 0,
                           ind[0], ind[1], ind[2], ind[3], ii);
      }

    }
  } else if (tab->num_param == 6) {

    if (is_xill_model(param->model_type)) {
      for (ii = 0; ii < spec->n_incl; ii++) {
        interp_6d_tab(tab, spec->flu[ii], spec->n_ener,
                      fac, nfac, ind);
      }
    } else {
      for (ii = 0; ii < spec->n_incl; ii++) {
        interp_6d_tab_incl(tab, spec->flu[ii], spec->n_ener,
                           fac, nfac, ind, ii);
      }
    }

  }

  free(fac);

  return spec;
}

char *getXilltableNameUsingAlternativeIfNotExisting(char *stdname, char *altname, int *status) {
  char *nthcompFileName = stdname;

  if (checkIfTableExists(nthcompFileName, status) == 0) {
    nthcompFileName = altname;
    if (is_debug_run()) {
      printf(" *** warning: did not find standard table %s, trying to use alternatively %s \n", stdname, altname);
    }
  }
  return nthcompFileName;
}
/** load the xillver table and return its filename **/
char *get_init_xillver_table(xillTable **tab, xillParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  if (is_dens_model(param->model_type) && (param->prim_type == PRIM_SPEC_ECUT)) {
    if (cached_xill_tab_dens == NULL) {
      init_xillver_table(XILLTABLE_DENS_FILENAME, &cached_xill_tab_dens, param, status);
    }
    *tab = cached_xill_tab_dens;
    return XILLTABLE_DENS_FILENAME;

  } else if (is_ns_model(param->model_type)) {
    if (cached_xill_tab_ns == NULL) {
      init_xillver_table(XILLTABLE_NS_FILENAME, &cached_xill_tab_ns, param, status);
      CHECK_STATUS_RET(*status, NULL);
    }
    *tab = cached_xill_tab_ns;
    return XILLTABLE_NS_FILENAME;

  } else if (is_co_model(param->model_type)) {
    if (cached_xill_tab_co == NULL) {
      init_xillver_table(XILLTABLE_CO_FILENAME, &cached_xill_tab_co, param, status);
      CHECK_STATUS_RET(*status, NULL);
    }
    *tab = cached_xill_tab_co;
    return XILLTABLE_CO_FILENAME;

  } else if (is_CpD_model(param->model_type)) {

    if (cached_xill_tab_dens_nthcomp == NULL) {
      init_xillver_table(XILLTABLE_NTHCOMP_FILENAME, &cached_xill_tab_dens_nthcomp, param, status);
      CHECK_STATUS_RET(*status, NULL);
    }
    *tab = cached_xill_tab_dens_nthcomp;
    return XILLTABLE_NTHCOMP_FILENAME;

  } else if (param->prim_type == PRIM_SPEC_NTHCOMP) {

    // currently we have a new and also an older table for the nthcomp version (one includes more
    // densities, but is not in the standard distribution
    char *nthcompFileName = getXilltableNameUsingAlternativeIfNotExisting(
        XILLTABLE_NTHCOMP_FILENAME, XILLTABLE_NTHCOMP_FILENAME_OLD, status);

    if (cached_xill_tab_nthcomp == NULL) {
      init_xillver_table(nthcompFileName, &cached_xill_tab_nthcomp, param, status);
    }

    *tab = cached_xill_tab_nthcomp;
    return nthcompFileName;

  } else {

    if (cached_xill_tab == NULL) {
      init_xillver_table(XILLTABLE_FILENAME, &cached_xill_tab, param, status);
    }
    *tab = cached_xill_tab;
    return XILLTABLE_FILENAME;
  }

}

/** the main routine for the xillver table: returns a spectrum for the given parameters
 *  - decides if the table needs to be initialized and/or more data loaded
 *  - automatically normalizes  the spectra to logN=1e15 and logXi=0
 * */
xill_spec *get_xillver_spectra(xillParam *param, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  xillTable *tab = NULL;
  char *fname = get_init_xillver_table(&tab, param, status);

  CHECK_STATUS_RET(*status, NULL);
  assert(fname != NULL);

  // =1=  get the inidices
  int *ind = xillInd_from_parInput(param, tab, status);

  // =2=  check if the necessary spectra are loaded (we only open the file once)
  check_xillTable_cache(fname, param, tab, ind, status);

  // =3= interpolate values
  xill_spec *spec = interp_xill_table(tab, param, ind, status);

  CHECK_RELXILL_DEFAULT_ERROR(status);

  free(ind);
  return spec;
}

void norm_xillver_spec(xill_spec *spec, double incl) {

  /** adds the proper flux normalization for a semi-infinate slab
   *  under inclination angle incl */
  int ii;
  for (ii = 0; ii < spec->n_ener; ii++) {
    spec->flu[0][ii] *= 0.5 * cos(incl * M_PI / 180);
  }

}

void free_xillTable(xillTable *tab) {
  if (tab != NULL) {

    int ii;
    if (tab->data_storage != NULL) {
      for (ii = 0; ii < tab->num_elements; ii++) {
        if (tab->data_storage[ii] != NULL) {
          free(tab->data_storage[ii]);
        }

      }
      free(tab->data_storage);
    }

    if (tab->param_vals != NULL) {
      for (ii = 0; ii < tab->num_param; ii++) {
        free(tab->param_vals[ii]);

        free(tab->param_names[ii]);
      }

    }
    free(tab->param_vals);
    free(tab->param_names);

    free(tab->num_param_vals);
    free(tab->param_index);

    free(tab->elo);
    free(tab->ehi);

    free(tab);
  }
}

void free_cached_xillTable(void) {
  free_xillTable(cached_xill_tab);
  free_xillTable(cached_xill_tab_dens);
  free_xillTable(cached_xill_tab_nthcomp);
}
