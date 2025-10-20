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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "xilltable_reading.h"
#include "relutility.h"
#include "common.h"


// possible parameters for the xillver tables
int global_param_index[] = {
  PARAM_GAM, PARAM_AFE, PARAM_LXI, PARAM_ECT, PARAM_KTE, PARAM_DNS,
  PARAM_KTB, PARAM_ACO, PARAM_FRA, PARAM_INC, PARAM_RID
};
char* global_param_names[] = {
  NAME_GAM, NAME_AFE, NAME_LXI, NAME_ECT, NAME_KTE, NAME_DNS,
  NAME_KTB, NAME_ACO, NAME_FRA, NAME_INC, NAME_RID
};

// storage for the tables
xillTable* cached_xill_tab = NULL;
xillTable* cached_xill_tab_dens = NULL;
xillTable* cached_xill_tab_nthcomp = NULL;
xillTable* cached_xill_tab_dens_nthcomp = NULL;
xillTable* cached_xill_tab_ns = NULL;
xillTable* cached_xill_tab_co = NULL;
xillTable* cached_xill_tab_rrad_bb = NULL;

static int get_num_elem(const int* n_parvals, int npar)
{
  assert(npar >= 1);
  long num_elem = n_parvals[0];

  for (int ii = 1; ii < npar; ii++)
  {
    num_elem *= n_parvals[ii];
  }

  assert(num_elem <= INT_MAX);

  return (int)num_elem;
}

static void init_xilltable_data_struct(xillTable* tab, int* status)
{
  CHECK_STATUS_VOID(*status);
  assert(tab != NULL);
  assert(tab->param_vals != NULL);
  assert(tab->data_storage == NULL);

  tab->num_elements = get_num_elem(tab->num_param_vals, tab->num_param);

  tab->data_storage = (float**)malloc(sizeof(float*) * tab->num_elements);
  CHECK_MALLOC_VOID_STATUS(tab->data_storage, status)

  // important to make sure everything is set to NULL (used to only load spectra if !=NULL)
  int ii;
  for (ii = 0; ii < tab->num_elements; ii++)
  {
    tab->data_storage[ii] = NULL;
  }
}

/** get a new and empty rel table (structure will be allocated)  */
xillTable* new_xillTable(int num_param, int* status)
{
  xillTable* tab = (xillTable*)malloc(sizeof(xillTable));
  CHECK_MALLOC_RET_STATUS(tab, status, tab)

  tab->num_param = num_param;


  tab->num_param_vals = (int*)malloc(sizeof(int) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->num_param_vals, status, NULL)

  tab->param_index = (int*)malloc(sizeof(int) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->param_index, status, NULL)

  tab->param_names = (char**)malloc(sizeof(char*) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->param_index, status, NULL)

  tab->param_vals = (float**)malloc(sizeof(float*) * num_param);
  CHECK_MALLOC_RET_STATUS(tab->param_vals, status, NULL)

  int ii;
  for (ii = 0; ii < num_param; ii++)
  {
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

int is_6dim_table(int model_type)
{
  if (is_co_model(model_type))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

static int get_num_param_auto(fitsfile* fptr, int* status)
{
  CHECK_STATUS_RET(*status, 0);

  long n;
  fits_get_num_rows(fptr, &n, status);
  relxill_check_fits_error(status);

  return (int)n;
}

/* routine to get the pointer to the pre-defined parameter index array
 * (i.e., with this we know which parameter relates to the input parameters) */
static void set_parindex_from_parname(int* pindex, char** pname, int n, int* status)
{
  int ii;
  // loop over given parameters
  for (ii = 0; ii < n; ii++)
  {
    pindex[ii] = -1;
    assert(pname!=NULL);

    // now loop over all possible parameters
    int jj;
    for (jj = 0; jj < N_PARAM_MAX; jj++)
    {
      if (strcmp(pname[ii], global_param_names[jj]) == 0)
      {
        pindex[ii] = global_param_index[jj];
        continue;
      }
    }

    if (pindex[ii] == -1)
    {
      RELXILL_ERROR(" parameter given in xillver table not found\n", status);
      printf("    parameter ** %s ** from xillver table, not known to relxill \n", pname[ii]);
      printf("    please make sure you downloaded the correct table \n");
    }
  }
}

void print_xilltable_parameters(const xillTable* tab, char* const * xilltab_parname)
{
  assert(tab != NULL);

  int ii;
  for (ii = 0; ii < tab->num_param; ii++)
  {
    printf(" loaded parameter %s  (index=%i) \t -  %02i values from %.2f to %.2f\n",
           xilltab_parname[ii], tab->param_index[ii], tab->num_param_vals[ii],
           tab->param_vals[ii][0], tab->param_vals[ii][tab->num_param_vals[ii] - 1]);
  }
}

// read the parameters of the xillver FITS table
static void get_xilltable_parameters(fitsfile* fptr, xillTable* tab, int* status)
{
  CHECK_STATUS_VOID(*status);

  int extver = 0;
  fits_movnam_hdu(fptr, BINARY_TBL, "PARAMETERS", extver, status);
  if (*status != EXIT_SUCCESS)
  {
    printf(" *** error moving to extension PARAMETERS in the xillver table\n");
    return;
  }

  // we know the column numbers
  int colnum_n = 9;

  long n;
  if (fits_get_num_rows(fptr, &n, status)) return;

  assert(tab != NULL);

  if (tab->num_param != n)
  {
    RELXILL_ERROR("wrong format of the xillver table (not enough or too many parameter values tabulated)", status);
  }

  int anynul = 0;
  double nullval = 0.0;

  int ii;
  char strnull[10];
  strcpy(strnull, " ");

  assert(tab->param_names != NULL);

  for (ii = 0; ii < tab->num_param; ii++)
  {
    tab->param_names[ii] = (char*)malloc(sizeof(char) * 16);
    CHECK_MALLOC_VOID_STATUS(tab->param_names[ii], status)
    memset(tab->param_names[ii], 0, 16); // Initialize to zero to ensure null termination
  }
  fits_read_col(fptr, TINT, colnum_n, 1, 1, tab->num_param, &nullval, tab->num_param_vals, &anynul, status);
  fits_read_col(fptr, TSTRING, 1, 1, 1, tab->num_param, strnull, tab->param_names, &anynul, status);

  // Ensure all strings are properly null-terminated (fits_read_col, might not guarantee this)
  for (ii = 0; ii < tab->num_param; ii++)
  {
    tab->param_names[ii][15] = '\0'; // Force null termination
  }

  for (ii = 0; ii < tab->num_param; ii++)
  {
    int colnum_vals = 10;
    /** the we load the parameter values **/
    float** ptr_val = &(tab->param_vals[ii]);

    if (ptr_val == NULL)
    {
      RELXILL_ERROR(" *** error loading the xillver parameter from the table", status);
      return;
    }
    (*ptr_val) = (float*)malloc(tab->num_param_vals[ii] * sizeof(float));
    CHECK_MALLOC_VOID_STATUS(*ptr_val, status)

    fits_read_col(fptr, TFLOAT, colnum_vals, ii + 1, 1, tab->num_param_vals[ii], &nullval, *ptr_val, &anynul,
                  status);

    CHECK_STATUS_BREAK(*status);
  }

  set_parindex_from_parname(tab->param_index, tab->param_names, tab->num_param, status);

  if (is_debug_run())
  {
    print_xilltable_parameters(tab, tab->param_names);
  }

  // now we set a pointer to the inclination separately (it is by definition the last parameter of each xillver table)
  assert(tab->param_names[tab->num_param - 1][0] == 'I'
    && tab->param_names[tab->num_param - 1][1] == 'n'
    && tab->param_names[tab->num_param - 1][2] == 'c'
    && tab->param_names[tab->num_param - 1][3] == 'l'
  );

  tab->incl = tab->param_vals[tab->num_param - 1];
  tab->n_incl = tab->num_param_vals[tab->num_param - 1];
}

static void get_xilltable_ener(int* n_ener, float** elo, float** ehi, fitsfile* fptr, int* status)
{
  CHECK_STATUS_VOID(*status);

  int extver = 0;
  fits_movnam_hdu(fptr, BINARY_TBL, "ENERGIES", extver, status);
  if (*status != EXIT_SUCCESS)
  {
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
  *n_ener = (int)nrows;

  *elo = (float*)malloc((*n_ener) * sizeof(float));
  CHECK_MALLOC_VOID_STATUS(*elo, status)
  *ehi = (float*)malloc((*n_ener) * sizeof(float));
  CHECK_MALLOC_VOID_STATUS(*ehi, status)

  int anynul = 0;
  double nullval = 0.0;
  LONGLONG nelem = (LONGLONG)(*n_ener);
  fits_read_col(fptr, TFLOAT, colnum_elo, 1, 1, nelem, &nullval, *elo, &anynul, status);
  fits_read_col(fptr, TFLOAT, colnum_ehi, 1, 1, nelem, &nullval, *ehi, &anynul, status);

  relxill_check_fits_error(status);

  CHECK_RELXILL_ERROR("reading of energy grid of the xillver table failed", status);
}

/**
 * @brief Convert xillver table parameters to a global parameter values array
 * @details Creates an array indexed by global parameter identifiers (PARAM_GAM, PARAM_AFE, etc.)
 *          containing the corresponding values from the xillTableParam structure. This allows
 *          unified access to all xillver parameters using their global index constants.
 * @param param Input structure containing xillver table parameters
 * @param status Status variable for error handling
 * @return Pointer to the allocated float array with N_PARAM_MAX elements, or NULL on error.
 *         Caller is responsible for freeing the returned array.
 */
float* get_xilltab_global_paramvals(const xillTableParam* param, int* status)
{
  CHECK_STATUS_RET(*status, NULL);

  float* global_param_vals = (float*)malloc(N_PARAM_MAX * sizeof(float));
  CHECK_MALLOC_RET_STATUS(global_param_vals, status, NULL)

  global_param_vals[PARAM_GAM] = (float)param->gam;
  global_param_vals[PARAM_AFE] = (float)param->afe;
  global_param_vals[PARAM_LXI] = (float)param->lxi;
  global_param_vals[PARAM_ECT] = (float)param->ect;
  global_param_vals[PARAM_DNS] = (float)param->dens;
  global_param_vals[PARAM_KTB] = (float)param->kTbb;
  global_param_vals[PARAM_FRA] = (float)param->frac_pl_bb;
  global_param_vals[PARAM_RID] = (float)param->radindex;
  global_param_vals[PARAM_INC] = (float)param->incl;

  return global_param_vals;
}

static int get_xillspec_rownum(const int* num_param_vals, int num_param, int nn, int ii, int jj, int kk, int ll, int mm)
{
  switch (num_param)
  {
  case 2:
    return (ii * num_param_vals[1] + jj + 1);

  case 5:
    return (((ii * num_param_vals[1] + jj) * num_param_vals[2]
      + kk) * num_param_vals[3] + ll) * num_param_vals[4] + mm + 1;

  case 6:
    return ((((nn * num_param_vals[1] + ii) * num_param_vals[2]
      + jj) * num_param_vals[3] + kk) * num_param_vals[4] + ll) * num_param_vals[5] + mm + 1;

  default:
    return -1;
  }
}

static void set_dat(float* spec, const xillTable* tab, int i0, int i1, int i2, int i3, int i4, int i5)
{
  int row_index = get_xillspec_rownum(tab->num_param_vals, tab->num_param,
                                      i0, i1, i2, i3, i4, i5);
  tab->data_storage[row_index - 1] = spec;
}

// get one Spectrum from the Data Storage
float* get_xillspec_table(const xillTable* tab, int i0, int i1, int i2, int i3, int i4, int i5)
{
  int rowindex = get_xillspec_rownum(tab->num_param_vals, tab->num_param,
                                     i0, i1, i2, i3, i4, i5);

  return tab->data_storage[rowindex - 1];
}

static char* getFullPathTableName(const char* filename, int* status)
{
  int MAXSIZE = 1000;
  char* fullfilename = (char*)malloc(sizeof(char) * MAXSIZE);
  CHECK_MALLOC_RET_STATUS(fullfilename, status, NULL)

  if (sprintf(fullfilename, "%s/%s", get_relxill_table_path(), filename) == -1)
  {
    RELXILL_ERROR("failed to construct full path the rel table", status);
    return NULL;
  }

  return fullfilename;
}

int checkIfTableExists(const char* filename, int* status)
{
  char* fullfilename = getFullPathTableName(filename, status);
  CHECK_STATUS_RET(*status, 0);

  fitsfile* fptr = NULL;
  int statusTableExists = EXIT_SUCCESS;
  int tableExists = 0;

  if (fits_open_table(&fptr, fullfilename, READONLY, &statusTableExists) == 0)
  {
    tableExists = 1;
  }

  if (fptr != NULL)
  {
    fits_close_file(fptr, &statusTableExists);
  }

  free(fullfilename);

  return tableExists;
}

fitsfile* open_fits_table_stdpath(const char* filename, int* status)
{
  CHECK_STATUS_RET(*status, NULL);

  char* full_filename = getFullPathTableName(filename, status);
  CHECK_STATUS_RET(*status, NULL);

  fitsfile* fptr = NULL;
  if (fits_open_table(&fptr, full_filename, READONLY, status))
  {
    RELXILL_ERROR("opening of the table failed", status);
    printf("    either the full path given (%s) is wrong \n", full_filename);
    printf("    or you need to download the table ** %s **  from \n", filename);
    printf("    https://www.sternwarte.uni-erlangen.de/research/relxill/ \n");
  }

  free(full_filename);

  return fptr;
}

void init_xillver_table(const char* filename, xillTable** inp_tab, int* status)
{
  CHECK_STATUS_VOID(*status);

  assert(*inp_tab == NULL); // only write the table if it is empty
  xillTable* tab = NULL;
  fitsfile* fptr = NULL;

  print_version_number();

  do
  {
    // loop to ensure if *status!=EXIT_SUCCESS no further functions are called
    fptr = open_fits_table_stdpath(filename, status);
    CHECK_STATUS_BREAK(*status);

    int num_param = get_num_param_auto(fptr, status);
    CHECK_STATUS_BREAK(*status);

    tab = new_xillTable(num_param, status);
    CHECK_STATUS_BREAK(*status);
    assert(tab != NULL);

    get_xilltable_ener(&(tab->n_ener), &(tab->elo), &(tab->ehi), fptr, status);
    CHECK_STATUS_BREAK(*status);

    get_xilltable_parameters(fptr, tab, status);
    CHECK_STATUS_BREAK(*status);

    init_xilltable_data_struct(tab, status);
    CHECK_STATUS_BREAK(*status);
  }
  while (0);

  if (*status == EXIT_SUCCESS && tab != NULL)
  {
    // assign the value
    *inp_tab = tab;
  }
  else
  {
    printf(" *** error *** initializing of the XILLVER table %s failed \n", filename);
    if (tab != NULL)
    {
      free_xillTable(&tab);
    }
  }

  if (fptr != NULL)
  {
    fits_close_file(fptr, status);
  }
}

/**
 * @brief: check if the given parameter index is in the xillver table
 * @return
 *  -1: not found
 *  ii: return index where this value is found in the param_index array
 **/
int get_xilltab_param_index(const xillTable* tab, int unique_param_identifier)
{
  int ii;
  for (ii = 0; ii < tab->num_param; ii++)
  {
    if (tab->param_index[ii] == unique_param_identifier)
    {
      return ii;
    }
  }
  return -1;
}

/* convert the index from the (actual) xillver FITS table parameter
/ extension to the code internal param_index structure, which leaves
/ the first bins empty to fill up the array until the last bin
*/
static int convertTo6dimTableIndex(int num_param, int indexTable)
{
  const int num_param_max = 6;
  assert(num_param <= 6 || num_param >= 2);
  if (num_param != num_param_max)
  {
    indexTable += num_param_max - num_param; // for 5dim, the first index is left empty
  }
  return indexTable;
}


void renorm_xill_spec(float* spec, int n, double lxi, double dens)
{
  for (int ii = 0; ii < n; ii++)
  {
    spec[ii] /= pow(10, lxi); // do not cast to float (fails refdata) NOLINT(*-narrowing-conversions)
    if (fabs(dens - 15) > 1e-6)
    {
      spec[ii] /= pow(10, dens - 15); // do not cast to float (fails refdata) NOLINT(*-narrowing-conversions)
    }
  }
}


static void normalizeXillverSpecLogxiDensity(float* spec,
                                             const xillTable* tab,
                                             double density,
                                             double logxi,
                                             const int* indVals)
{
  int ind_dens = get_xilltab_param_index(tab, PARAM_DNS);
  int ind_lxi = get_xilltab_param_index(tab, PARAM_LXI);


  // ind==-1 if parameter not contained
  if (ind_dens >= 0)
  {
    int ind_dens_6dim = convertTo6dimTableIndex(tab->num_param, ind_dens);
    density = tab->param_vals[ind_dens][indVals[ind_dens_6dim]];
  }
  if (ind_lxi >= 0)
  {
    int ind_lxi_6dim = convertTo6dimTableIndex(tab->num_param, ind_lxi);
    logxi = tab->param_vals[ind_lxi][indVals[ind_lxi_6dim]];
  }

  renorm_xill_spec(spec, tab->n_ener, logxi, density);
}


void xilltable_fits_load_single_spec(const char* fname,
                                     fitsfile** fptr,
                                     const xillTable* tab,
                                     double defDensity,
                                     double defLogxi,
                                     int nn,
                                     int ii,
                                     int jj,
                                     int kk,
                                     int ll,
                                     int mm,
                                     int* status)
{
  CHECK_STATUS_VOID(*status);

  // open the FITS file if not already open
  if (*fptr == NULL)
  {
    *fptr = open_fits_table_stdpath(fname, status);
    CHECK_STATUS_VOID(*status);
  }

  int extver = 0;
  fits_movnam_hdu(*fptr, BINARY_TBL, "SPECTRA", extver, status);
  if (*status != EXIT_SUCCESS)
  {
    printf(" *** error moving to extension SPECTRA in the xillver table\n");
    return;
  }

  float* spec = (float*)malloc(tab->n_ener * sizeof(float));
  CHECK_MALLOC_VOID_STATUS(spec, status)

  int colnum_spec = 2;
  int anynul = 0;
  double nullval = 0.0;
  LONGLONG nelem = (LONGLONG)tab->n_ener;

  // get the row number (this is how the xillver table is defined)
  int rownum = get_xillspec_rownum(tab->num_param_vals, tab->num_param, nn, ii, jj, kk, ll, mm);

  fits_read_col(*fptr, TFLOAT, colnum_spec, rownum, 1, nelem, &nullval, spec, &anynul, status);
  if (*status != EXIT_SUCCESS)
  {
    printf("\n *** ERROR *** failed reading table %s  (rownum %i) \n",
           fname, rownum);
    relxill_check_fits_error(status);
    return;
  }

  int indexArray[] = {nn, ii, jj, kk, ll, mm};
  normalizeXillverSpecLogxiDensity(spec, tab, defDensity, defLogxi, indexArray);

  set_dat(spec, tab, nn, ii, jj, kk, ll, mm);
}


enum xillTableIds get_xilltable_id(int model_id, int prim_type)
{
  if (is_ns_model(model_id))
  {
    return XILLTABLE_ID_NS;
  }

  if (is_co_model(model_id))
  {
    return XILLTABLE_ID_CO;
  }

  if (prim_type == PRIM_SPEC_NTHCOMP)
  {
    return XILLTABLE_ID_NTHCOMP;
  }

  if (is_returnrad_bb_model(model_id))
  {
    if (model_id == MOD_TYPE_RELXILLBBRET)
    {
      // relxillBB uses xillverNS table
      return XILLTABLE_ID_NS;
    }
    return XILLTABLE_ID_RRAD; // relxillBBxill uses xillverRTR table
  }

  return XILLTABLE_ID_STANDARD;
}

xillTableParam* get_xilltab_param(const xillParam* param, int* status)
{
  CHECK_STATUS_RET(*status, NULL);

  xillTableParam* tab_param = malloc(sizeof(xillTableParam));
  CHECK_MALLOC_RET_STATUS(tab_param, status, NULL)

  tab_param->gam = param->gam;
  tab_param->afe = param->afe;
  tab_param->lxi = param->lxi;
  tab_param->ect = param->ect;
  tab_param->dens = param->dens;
  tab_param->kTbb = param->kTbb;
  tab_param->frac_pl_bb = param->frac_pl_bb;
  tab_param->incl = param->incl;
  tab_param->radindex = param->radzone;

  tab_param->model_type = param->model_type;
  tab_param->prim_type = param->prim_type;

  return tab_param;
}

/** load the xillver table and return its filename **/
const char* get_init_xillver_table(xillTable** tab, int model_type, int prim_type, int* status)
{
  CHECK_STATUS_RET(*status, NULL);

  switch (get_xilltable_id(model_type, prim_type))
  {
  case XILLTABLE_ID_STANDARD:
    {
      if (cached_xill_tab == NULL)
      {
        init_xillver_table(XILLTABLE_FILENAME, &cached_xill_tab, status);
      }
      *tab = cached_xill_tab;
      return XILLTABLE_FILENAME;
    }

  case XILLTABLE_ID_NTHCOMP:
    {
      if (cached_xill_tab_dens_nthcomp == NULL)
      {
        init_xillver_table(XILLTABLE_NTHCOMP_FILENAME, &cached_xill_tab_dens_nthcomp, status);
        CHECK_STATUS_RET(*status, NULL);
      }
      *tab = cached_xill_tab_dens_nthcomp;
      return XILLTABLE_NTHCOMP_FILENAME;
    }

  case XILLTABLE_ID_NS:
    {
      if (cached_xill_tab_ns == NULL)
      {
        init_xillver_table(XILLTABLE_NS_FILENAME, &cached_xill_tab_ns, status);
        CHECK_STATUS_RET(*status, NULL);
      }
      *tab = cached_xill_tab_ns;
      return XILLTABLE_NS_FILENAME;
    }

  case XILLTABLE_ID_CO:
    {
      if (cached_xill_tab_co == NULL)
      {
        // not working anymore as param does not exist
        //        assert(param->dens == 17); // the CO_Table is explicitly calculated for logN=17
        //        assert(param->lxi == 0); // the CO_Table does not have ionization as a parameter (weak test to assert it has its default value)
        init_xillver_table(XILLTABLE_CO_FILENAME, &cached_xill_tab_co, status);
        CHECK_STATUS_RET(*status, NULL);
      }
      *tab = cached_xill_tab_co;
      return XILLTABLE_CO_FILENAME;
    }

  case XILLTABLE_ID_RRAD:
    {
      if (cached_xill_tab_rrad_bb == NULL)
      {
        init_xillver_table(XILLTABLE_RRAD_FILENAME, &cached_xill_tab_rrad_bb, status);
        CHECK_STATUS_RET(*status, NULL);
      }
      *tab = cached_xill_tab_rrad_bb;
      return XILLTABLE_RRAD_FILENAME;
    }

  default:
    {
      RELXILL_ERROR("could not identify a xillver table for this model", status);
      return NULL;
    }
  }
}


void free_xillTable(xillTable** ptr_tab)
{
  xillTable* tab = *ptr_tab;
  if (tab != NULL)
  {
    int ii;
    if (tab->data_storage != NULL)
    {
      for (ii = 0; ii < tab->num_elements; ii++)
      {
        if (tab->data_storage[ii] != NULL)
        {
          free(tab->data_storage[ii]);
        }
      }
      free(tab->data_storage);
    }

    if (tab->param_vals != NULL)
    {
      for (ii = 0; ii < tab->num_param; ii++)
      {
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
    tab = NULL;
  }
}

void free_all_cached_xillTables(void)
{
  if (cached_xill_tab != NULL)
  {
    free_xillTable(&cached_xill_tab);
  }
  if (cached_xill_tab_dens != NULL)
  {
    free_xillTable(&cached_xill_tab_dens);
  }
  if (cached_xill_tab_nthcomp != NULL)
  {
    free_xillTable(&cached_xill_tab_nthcomp);
  }
  if (cached_xill_tab_dens_nthcomp != NULL)
  {
    free_xillTable(&cached_xill_tab_dens_nthcomp);
  }
  if (cached_xill_tab_ns != NULL)
  {
    free_xillTable(&cached_xill_tab_ns);
  }
  if (cached_xill_tab_co != NULL)
  {
    free_xillTable(&cached_xill_tab_co);
  }
  if (cached_xill_tab_rrad_bb != NULL)
  {
    free_xillTable(&cached_xill_tab_rrad_bb);
  }
}
