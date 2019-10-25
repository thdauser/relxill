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

    Copyright 2019 Thomas Dauser, Remeis Observatory & ECAP
*/

#include <zconf.h>
#include "xilltable.h"

// has to be 6DIM, where are the parameters
// [param->gam, param->afe, param->lxi, param->ect, param->dens, param->incl]
int XILL_PARAM_INDEX[] = {0, 1, 2, 3, -1, 4};
int XILL_DENS_PARAM_INDEX[] = {0, 1, 2, -1, 3, 4};
int XILL_NTHCOMP_PARAM_INDEX[] = {0, 1, 2, -1, 3, 4};

xillTable *cached_xill_tab = NULL;
xillTable *cached_xill_tab_dens = NULL;
xillTable *cached_xill_tab_nthcomp = NULL;

static void init_5dim_xilltable(float ******dat, xillTable *tab, int *status) {
// allocate the full arrays just with pointers; make sure everything is set to NULL in the end

    int istart = (tab->num_param + 1) - MAX_DIM_TABLE; // if 5dim, we start at 0, for 6dim start at 1

    dat = (float ******) malloc(sizeof(float *****) * tab->num_param_vals[istart]);
    CHECK_MALLOC_VOID_STATUS(dat, status)

    int ii;
    int jj;
    int kk;
    int ll;
    int mm;
    for (ii = 0; ii < tab->num_param_vals[istart]; ii++) {
        dat[ii] = (float *****) malloc(sizeof(float ****) * tab->num_param_vals[istart + 1]);
        CHECK_MALLOC_VOID_STATUS(dat[ii], status)

        for (jj = 0; jj < tab->num_param_vals[istart + 1]; jj++) {
            dat[ii][jj] = (float ****) malloc(sizeof(float ***) * tab->num_param_vals[istart + 2]);
            CHECK_MALLOC_VOID_STATUS(dat[ii][jj], status)

            for (kk = 0; kk < tab->num_param_vals[istart + 2]; kk++) {
                dat[ii][jj][kk] = (float ***) malloc(sizeof(float **) * tab->num_param_vals[istart + 3]);
                CHECK_MALLOC_VOID_STATUS(dat[ii][jj][kk], status)

                for (ll = 0; ll < tab->num_param_vals[istart + 3]; ll++) {
                    dat[ii][jj][kk][ll] = (float **) malloc(sizeof(float *) * tab->num_param_vals[istart + 4]);
                    CHECK_MALLOC_VOID_STATUS(dat[ii][jj][kk][ll], status)

                    for (mm = 0; mm < tab->num_param_vals[istart + 4]; mm++) {
                        dat[ii][jj][kk][ll][mm] = NULL;
                    }
                }
            }
        }
    }
}


static void init_xilltable_data_struct(xillTable *tab, int *status) {

    CHECK_STATUS_VOID(*status)

    // check if we have only a 5dim table, then we can set this to make the rest of the code easier
    int npar0 = tab->num_param_vals[0];
    if (tab->num_param == 5) {
        npar0 = 1;
    }

    tab->data_storage = (float *******) malloc(sizeof(float ******) * npar0);
    CHECK_MALLOC_VOID_STATUS(tab->data_storage, status)

    int ii;
    for (ii = 0; ii < npar0; ii++) {
        init_5dim_xilltable(tab->data_storage[ii], tab, status);
    }


    if (tab->num_param == 5) {
        tab->dat = (tab->data_storage[0]);  // set the pointer to the 5D table
    } else {
        tab->dat = (tab->data_storage);    // set the pointer to the 6D table
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

    tab->param_vals = (float **) malloc(sizeof(float *) * num_param);
    CHECK_MALLOC_RET_STATUS(tab->param_vals, status, NULL)

    int ii;
    for (ii = 0; ii < num_param; ii++) {
        tab->param_vals[ii] = NULL;
    }


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

static int is_6dim_table(int model_type) {
    if (model_type == MOD_TYPE_RELXILLCO) {
        return 1;
    } else {
        return 0;
    }
}

static int get_num_param(xillParam *param) {

    int model_type = param->model_type;

    int num_param = XILLTABLE_N_PARAM;
    if (is_6dim_table(model_type)) {
        num_param = 6;
    }

    return num_param;
}

/* routine to get the pointer to the pre-defined parameter index array
 * (i.e. with this we know which parameter relates to the input parameters )
 * TODO: in principle this routine is not necessary and the information can be taken from the table */
static int *get_ptr_param_index(xillParam *param) {
    int *param_index;

    if (is_dens_model(param->model_type)) {
        param_index = XILL_DENS_PARAM_INDEX;
    } else if (param->prim_type == PRIM_SPEC_NTHCOMP) {
        param_index = XILL_NTHCOMP_PARAM_INDEX;
    } else {
        param_index = XILL_PARAM_INDEX;
    }
    return param_index;
}


/** read the parameters of the xillver FITS table   */
static void get_xilltable_parameters(fitsfile *fptr, xillTable *tab, xillParam *param, int *status) {

    CHECK_STATUS_VOID(*status)

    int extver = 0;
    fits_movnam_hdu(fptr, BINARY_TBL, "PARAMETERS", extver, status);
    if (*status != EXIT_SUCCESS) {
        printf(" *** error moving to extension PARAMETERS in the xillver table\n");
        return;
    }

    // first set the number of parameters to the table
    tab->num_param = get_num_param(param);

    // get a pointer to the parameter index array
    tab->param_index = get_ptr_param_index(param);

    // initialize memory to store the number of values of each parameter
    tab->num_param_vals = (int *) malloc(sizeof(int) * tab->num_param);
    CHECK_MALLOC_VOID_STATUS(tab->num_param_vals, status)

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
    /* allocate space for string column value */
    char **xilltab_parname = (char **) malloc(sizeof(char *) * tab->num_param);
    CHECK_MALLOC_VOID_STATUS(xilltab_parname, status)


    strcpy(strnull, " ");

    /** also get the name **/
    // TODO: get reading a string working
    for (ii = 0; ii < tab->num_param; ii++) {
        xilltab_parname[ii] = (char *) malloc(sizeof(char) * 32);
        CHECK_MALLOC_VOID_STATUS(xilltab_parname[ii], status)
    }
    fits_read_col(fptr, TINT, colnum_n, 1, 1, tab->num_param, &nullval, tab->num_param_vals, &anynul, status);
    fits_read_col(fptr, TSTRING, 1, 1, 1, tab->num_param, strnull, xilltab_parname, &anynul, status);


    for (ii = 0; ii < tab->num_param; ii++) {
        /** first get the number of parameters and check if those are correct **/
        //   fits_read_col(fptr, TINT, colnum_n, ii + 1, 1, 1, &nullval, &(tab->num_param_vals[ii]), &anynul, status);

        //      tab->param_vals[ii] = (float*) malloc(sizeof(float) * tab->num_param_vals[ii]);

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


        if (is_debug_run()) {
            printf(" loading parameter %s  \t -  %02i values from %.2f to %.2f\n",
                   xilltab_parname[ii], tab->num_param_vals[ii],
                   tab->param_vals[ii][0], tab->param_vals[ii][tab->num_param_vals[ii] - 1]);
        }


        CHECK_STATUS_BREAK(*status)
    }

    // now we set a pointer to the inclination separately (assuming it is the last parameter)
    tab->incl = tab->param_vals[tab->num_param - 1];
    tab->n_incl = tab->num_param_vals[tab->num_param - 1];

}



// static void get_reltable_axis(int nrows, float** val, char* extname, char* colname, fitsfile* fptr, int* status){
static void get_xilltable_ener(int *n_ener, float **elo, float **ehi, fitsfile *fptr, int *status) {

    CHECK_STATUS_VOID(*status)

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
    LONGLONG nelem = (LONGLONG)(*n_ener);
    fits_read_col(fptr, TFLOAT, colnum_elo, 1, 1, nelem, &nullval, *elo, &anynul, status);
    fits_read_col(fptr, TFLOAT, colnum_ehi, 1, 1, nelem, &nullval, *ehi, &anynul, status);

    relxill_check_fits_error(status);

    // simply let's check if something went wrong
    CHECK_RELXILL_ERROR("reading of energy grid of the xillver table failed", status);

}

static float *get_xill_param_vals_array(xillParam *param, int *status) {

    CHECK_STATUS_RET(*status, NULL)

    float *param_vals = (float *) malloc(N_PARAM_MAX * sizeof(float));
    CHECK_MALLOC_RET_STATUS(param_vals, status, NULL)
    // store the input in a variable
    param_vals[PARAM_GAM] = (float) param->gam;
    param_vals[PARAM_AFE] = (float) param->afe;
    param_vals[PARAM_LXI] = (float) param->lxi;
    param_vals[PARAM_ECT] = (float) param->ect;
    param_vals[PARAM_DNS] = (float) param->dens;
    param_vals[PARAM_INC] = (float) param->incl;

    return param_vals;
}

static int *get_xill_indices(xillParam *param, xillTable *tab, int *status) {

    CHECK_STATUS_RET(*status, NULL)
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
    return ind;
}


static fitsfile *open_xillver_tab(char *filename, int *status) {

    CHECK_STATUS_RET(*status, NULL)

    fitsfile *fptr = NULL;
    char fullfilename[500];
    // get the full filename
    if (sprintf(fullfilename, "%s/%s", get_relxill_table_path(), filename) == -1) {
        RELXILL_ERROR("failed to construct full path the rel table", status);
        return NULL;
    }

    // open the file
    if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
        CHECK_RELXILL_ERROR("opening of the xillver table failed", status);
        printf("    either the full path given (%s) is wrong \n", fullfilename);
        printf("    or you need to download the table ** %s **  from \n", filename);
        printf("    http://www.sternwarte.uni-erlangen.de/research/relxill/ \n");
        return NULL;
    }

    // free(fullfilename);

    return fptr;
}

/** load the complete relline table */
void init_xillver_table(char *filename, xillTable **inp_tab, xillParam *param, int *status) {

    CHECK_STATUS_VOID(*status)

    xillTable *tab = (*inp_tab);
    fitsfile *fptr = NULL;

    print_version_number(status);

    fptr = open_xillver_tab(filename, status);

    assert(tab == NULL);

    int num_param = get_num_param(param);  // need a routine to get the number of parameters here

    /** =1= allocate space for the new table  **/
    tab = new_xillTable(num_param, status);

    /** =2= now load the energy grid **/
    get_xilltable_ener(&(tab->n_ener), &(tab->elo), &(tab->ehi), fptr, status);

    /** =3= and now the stored parameter values (also check if the correct number of parameters) **/
    get_xilltable_parameters(fptr, tab, param, status);

    /** =4= finally set up the complete data structure **/
    init_xilltable_data_struct(tab, status);

    if (*status == EXIT_SUCCESS) {
        // should be set by previous routine
        assert(tab != NULL);
        assert(tab->elo != NULL);
        assert(tab->data_storage != NULL);
        // assigne the value
        (*inp_tab) = tab;
    } else {
        free_xillTable(tab);
    }

    if (fptr != NULL) { fits_close_file(fptr, status); }

}

static int
get_xillspec_rownum(const int *num_param_vals, int num_param, int nn, int ii, int jj, int kk, int ll, int mm) {

    if (num_param == 5) {
        return (((ii * num_param_vals[1] + jj) * num_param_vals[2]
                 + kk) * num_param_vals[3] + ll) * num_param_vals[4] + mm + 1;
    } else if (num_param == 6) {
        return ((((nn * num_param_vals[1] + ii) * num_param_vals[2]
                  + jj) * num_param_vals[3] + kk) * num_param_vals[4] + ll) * num_param_vals[5] + mm + 1;
        return -1;
    }

    return -1;
}


static void
load_single_spec(char *fname, fitsfile **fptr, xillTable *tab, int nn, int ii, int jj, int kk, int ll, int mm,
                 int *status) {

    CHECK_STATUS_VOID(*status)

    // open the fits file if not already open
    if (*fptr == NULL) {
        *fptr = open_xillver_tab(fname, status);
        CHECK_STATUS_VOID(*status)
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
    CHECK_STATUS_VOID(*status)

    if (tab->num_param == 5) {
        tab->data_storage[0][ii][jj][kk][ll][mm] = spec;
    } else {
        RELXILL_ERROR(" 6dim not implemented yet", status);
    }

}


void norm_xillver_spec(xill_spec *spec, double incl) {

    /** adds the proper flux normalization for a semi-infinate slab
     *  under inclination angle incl */
    int ii;
    for (ii = 0; ii < spec->n_ener; ii++) {
        spec->flu[0][ii] *= 0.5 * cos(incl * M_PI / 180);
    }

}

static void check_xillTable_cache(char *fname, xillTable *tab, const int *ind, int *status) {

    CHECK_STATUS_VOID(*status)

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
    // for 6dim
    if (tab->num_param == 6) {
        i0lo = ind[0];
        i0hi = ind[0] + 1;
    }


    for (nn = i0lo; nn <= i0hi; nn++) { // for 5dim this is a dummy loop
        for (ii = ind[istart]; ii <= ind[istart] + 1; ii++) {
            for (jj = ind[istart + 1]; jj <= ind[istart + 1] + 1; jj++) {
                for (kk = ind[istart + 2]; kk <= ind[istart + 2] + 1; kk++) {
                    for (ll = ind[istart + 3]; ll <= ind[istart + 3] + 1; ll++) {
                        // always load **all** incl bins as for relxill we will certainly need it
                        for (mm = 0; mm < tab->n_incl; mm++) {

                            if (tab->data_storage[nn][ii][jj][kk][ll][mm] == NULL) {
                                load_single_spec(fname, &fptr, tab, nn, ii, jj, kk, ll, mm, status);
                                CHECK_STATUS_VOID(*status)
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

    CHECK_STATUS_RET(*status, NULL)

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

static void interp_5d_tab_incl(float ******dat, double *flu, int n_ener,
                               double f1, double f2, double f3, double f4, int i1, int i2, int i3, int i4, int i5) {

    int ii;
    for (ii = 0; ii < n_ener; ii++) {
        flu[ii] =
                ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * dat[i1][i2][i3][i4][i5][ii] +
                 (f1) * (1.0 - f2) * (1.0 - f3) * dat[i1 + 1][i2][i3][i4][i5][ii] +
                 (1.0 - f1) * (f2) * (1.0 - f3) * dat[i1][i2 + 1][i3][i4][i5][ii] +
                 (1.0 - f1) * (1.0 - f2) * (f3) * dat[i1][i2][i3 + 1][i4][i5][ii] +
                 (f1) * (f2) * (1.0 - f3) * dat[i1 + 1][i2 + 1][i3][i4][i5][ii] +
                 (f1) * (1.0 - f2) * (f3) * dat[i1 + 1][i2][i3 + 1][i4][i5][ii] +
                 (1.0 - f1) * (f2) * (f3) * dat[i1][i2 + 1][i3 + 1][i4][i5][ii] +
                 (f1) * (f2) * (f3) * dat[i1 + 1][i2 + 1][i3 + 1][i4][i5][ii])
                * (1 - f4) +
                ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * dat[i1][i2][i3][i4 + 1][i5][ii] +
                 (f1) * (1.0 - f2) * (1.0 - f3) * dat[i1 + 1][i2][i3][i4 + 1][i5][ii] +
                 (1.0 - f1) * (f2) * (1.0 - f3) * dat[i1][i2 + 1][i3][i4 + 1][i5][ii] +
                 (1.0 - f1) * (1.0 - f2) * (f3) * dat[i1][i2][i3 + 1][i4 + 1][i5][ii] +
                 (f1) * (f2) * (1.0 - f3) * dat[i1 + 1][i2 + 1][i3][i4 + 1][i5][ii] +
                 (f1) * (1.0 - f2) * (f3) * dat[i1 + 1][i2][i3 + 1][i4 + 1][i5][ii] +
                 (1.0 - f1) * (f2) * (f3) * dat[i1][i2 + 1][i3 + 1][i4 + 1][i5][ii] +
                 (f1) * (f2) * (f3) * dat[i1 + 1][i2 + 1][i3 + 1][i4 + 1][i5][ii])
                * (f4);
    }

}

static void interp_5d_tab(float ******dat, double *flu, int n_ener,
                          double f1, double f2, double f3, double f4, double f5, int i1, int i2, int i3, int i4,
                          int i5) {

    int ii;
    for (ii = 0; ii < n_ener; ii++) {
        flu[ii] =
                (((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * dat[i1][i2][i3][i4][i5][ii] +
                  (f1) * (1.0 - f2) * (1.0 - f3) * dat[i1 + 1][i2][i3][i4][i5][ii] +
                  (1.0 - f1) * (f2) * (1.0 - f3) * dat[i1][i2 + 1][i3][i4][i5][ii] +
                  (1.0 - f1) * (1.0 - f2) * (f3) * dat[i1][i2][i3 + 1][i4][i5][ii] +
                  (f1) * (f2) * (1.0 - f3) * dat[i1 + 1][i2 + 1][i3][i4][i5][ii] +
                  (f1) * (1.0 - f2) * (f3) * dat[i1 + 1][i2][i3 + 1][i4][i5][ii] +
                  (1.0 - f1) * (f2) * (f3) * dat[i1][i2 + 1][i3 + 1][i4][i5][ii] +
                  (f1) * (f2) * (f3) * dat[i1 + 1][i2 + 1][i3 + 1][i4][i5][ii])
                 * (1 - f4) +
                 ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * dat[i1][i2][i3][i4 + 1][i5][ii] +
                  (f1) * (1.0 - f2) * (1.0 - f3) * dat[i1 + 1][i2][i3][i4 + 1][i5][ii] +
                  (1.0 - f1) * (f2) * (1.0 - f3) * dat[i1][i2 + 1][i3][i4 + 1][i5][ii] +
                  (1.0 - f1) * (1.0 - f2) * (f3) * dat[i1][i2][i3 + 1][i4 + 1][i5][ii] +
                  (f1) * (f2) * (1.0 - f3) * dat[i1 + 1][i2 + 1][i3][i4 + 1][i5][ii] +
                  (f1) * (1.0 - f2) * (f3) * dat[i1 + 1][i2][i3 + 1][i4 + 1][i5][ii] +
                  (1.0 - f1) * (f2) * (f3) * dat[i1][i2 + 1][i3 + 1][i4 + 1][i5][ii] +
                  (f1) * (f2) * (f3) * dat[i1 + 1][i2 + 1][i3 + 1][i4 + 1][i5][ii])
                 * (f4)) * (1.0 - f5) +
                (((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * dat[i1][i2][i3][i4][i5 + 1][ii] +
                  (f1) * (1.0 - f2) * (1.0 - f3) * dat[i1 + 1][i2][i3][i4][i5 + 1][ii] +
                  (1.0 - f1) * (f2) * (1.0 - f3) * dat[i1][i2 + 1][i3][i4][i5 + 1][ii] +
                  (1.0 - f1) * (1.0 - f2) * (f3) * dat[i1][i2][i3 + 1][i4][i5 + 1][ii] +
                  (f1) * (f2) * (1.0 - f3) * dat[i1 + 1][i2 + 1][i3][i4][i5 + 1][ii] +
                  (f1) * (1.0 - f2) * (f3) * dat[i1 + 1][i2][i3 + 1][i4][i5 + 1][ii] +
                  (1.0 - f1) * (f2) * (f3) * dat[i1][i2 + 1][i3 + 1][i4][i5 + 1][ii] +
                  (f1) * (f2) * (f3) * dat[i1 + 1][i2 + 1][i3 + 1][i4][i5 + 1][ii])
                 * (1 - f4) +
                 ((1.0 - f1) * (1.0 - f2) * (1.0 - f3) * dat[i1][i2][i3][i4 + 1][i5 + 1][ii] +
                  (f1) * (1.0 - f2) * (1.0 - f3) * dat[i1 + 1][i2][i3][i4 + 1][i5 + 1][ii] +
                  (1.0 - f1) * (f2) * (1.0 - f3) * dat[i1][i2 + 1][i3][i4 + 1][i5 + 1][ii] +
                  (1.0 - f1) * (1.0 - f2) * (f3) * dat[i1][i2][i3 + 1][i4 + 1][i5 + 1][ii] +
                  (f1) * (f2) * (1.0 - f3) * dat[i1 + 1][i2 + 1][i3][i4 + 1][i5 + 1][ii] +
                  (f1) * (1.0 - f2) * (f3) * dat[i1 + 1][i2][i3 + 1][i4 + 1][i5 + 1][ii] +
                  (1.0 - f1) * (f2) * (f3) * dat[i1][i2 + 1][i3 + 1][i4 + 1][i5 + 1][ii] +
                  (f1) * (f2) * (f3) * dat[i1 + 1][i2 + 1][i3 + 1][i4 + 1][i5 + 1][ii])
                 * (f4)) * (f5);

    }

}

/* cheap man's version of a 6dim interpolation */
static void interp_6d_tab_incl(float *******dat, double *flu, int n_ener,
                               double *fac, int nfac, int *ind, int iincl) {

    // we make sure this is used only for a 6dim table
    assert(nfac == 6);

    double s1[n_ener];
    double s2[n_ener];

    interp_5d_tab_incl(dat[ind[0]], s1, n_ener,
                       fac[1], fac[2], fac[3], fac[4],
                       ind[1], ind[2], ind[3], ind[4], iincl);

    interp_5d_tab_incl(dat[ind[0] + 1], s2, n_ener,
                       fac[1], fac[2], fac[3], fac[4],
                       ind[1], ind[2], ind[3], ind[4], iincl);

    int ii;
    for (ii = 0; ii < n_ener; ii++) {
        flu[ii] = interp_lin_1d(fac[0], s1[ii], s2[ii]);
    }
}

static xill_spec *interp_xill_table(xillTable *tab, xillParam *param, int *ind, int *status) {

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

    float *param_vals = get_xill_param_vals_array(param, status);

    int nfac = tab->num_param;
    double *fac = (double *) malloc(sizeof(double) * nfac);
    CHECK_MALLOC_RET_STATUS(fac, status, NULL)

    /* calculate interpolation factor for all parameters
     * ([nfac-1] is inclination, which might not be used ) */
    int pind;
    for (ii = 0; ii < nfac; ii++) {
        pind = tab->param_index[ii];
        fac[ii] = (param_vals[pind] - tab->param_vals[pind][ind[ii]]) /
                  (tab->param_vals[pind][ind[ii] + 1] - tab->param_vals[pind][ind[ii]]);
    }

    // can happen due to grav. redshift, although actually observed ecut is larger

    if (param->ect <= tab->param_vals[tab->param_index[PARAM_ECT]][0]) {
        fac[tab->param_index[PARAM_ECT]] = 0.0;
    }
    // can happen due to grav. redshift, although actually observed ecut is larger
    if (param->ect >=
        tab->param_vals[tab->param_index[PARAM_ECT]][tab->num_param_vals[tab->param_index[PARAM_ECT]] - 1]) {
        fac[tab->param_index[PARAM_ECT]] = 1.0;
    }

    if (tab->num_param == 5) {
        if (is_xill_model(param->model_type)) {
            assert(nfac == 5);
            interp_5d_tab(tab->dat, spec->flu[0], spec->n_ener,
                          fac[0], fac[1], fac[2], fac[3], fac[4],
                          ind[0], ind[1], ind[2], ind[3], ind[4]);


        } else {
            // do not interpolate over the inclination (so skip the last parameter)
            nfac--;
            assert(nfac == 4);
            // get the spectrum for EACH flux bin
            for (ii = 0; ii < spec->n_incl; ii++) {
                interp_5d_tab_incl(tab->dat, spec->flu[ii], spec->n_ener,
                                   fac[0], fac[1], fac[2], fac[3],
                                   ind[0], ind[1], ind[2], ind[3], ii);
            }

        }
    } else if (tab->num_param == 6) {

        if (is_xill_model(param->model_type)) {
            printf(" **** 6DIM Table not implemented for xillver model \n");
            return NULL;
        }

        assert(nfac == 6);
        for (ii = 0; ii < spec->n_incl; ii++) {
            interp_6d_tab_incl(tab->dat, spec->flu[ii], spec->n_ener,
                               fac, nfac, ind, ii);
        }

    }

    return spec;
}


/** load the xillver table and return its filename **/
char *get_init_xillver_table(xillTable **tab, xillParam *param, int *status) {

    CHECK_STATUS_RET(*status, NULL)

    if (is_dens_model(param->model_type)) {
        if (cached_xill_tab_dens == NULL) {
            init_xillver_table(XILLTABLE_DENS_FILENAME, &cached_xill_tab_dens, param, status);
        }
        *tab = cached_xill_tab_dens;
        return XILLTABLE_DENS_FILENAME;

    } else if (param->prim_type == PRIM_SPEC_NTHCOMP) {
        if (cached_xill_tab_nthcomp == NULL) {
            init_xillver_table(XILLTABLE_NTHCOMP_FILENAME, &cached_xill_tab_nthcomp, param, status);
        }
        *tab = cached_xill_tab_nthcomp;
        return XILLTABLE_NTHCOMP_FILENAME;

    } else {

        if (cached_xill_tab == NULL) {
            init_xillver_table(XILLTABLE_FILENAME, &cached_xill_tab, param, status);
        }
        *tab = cached_xill_tab;
        return XILLTABLE_FILENAME;
    }

}

/** the main routine for the xillver table: returns a spectrum for the given parameters
 *  (decides if the table needs to be initialized and/or more data loaded          */
xill_spec *get_xillver_spectra(xillParam *param, int *status) {

    CHECK_STATUS_RET(*status, NULL)

    xillTable *tab = NULL;
    char *fname = get_init_xillver_table(&tab, param, status);

    assert(tab != NULL);
    assert(fname != NULL);

    // =1=  get the inidices
    int *ind = get_xill_indices(param, tab, status);

    // =2=  check if the necessary spectra are loaded (we only open the file once)
    check_xillTable_cache(fname, tab, ind, status);

    // =3= interpolate values
    xill_spec *spec = interp_xill_table(tab, param, ind, status);

    free(ind);
    return spec;
}


static void free_xillTable_5dim(float ******dat, int n1, int n2, int n3, int n4, int n5) {

    int ii;
    int jj;
    int kk;
    int ll;
    int mm;

    if (dat != NULL) {
        for (ii = 0; ii < n1; ii++) {
            if (dat[ii] != NULL) {
                for (jj = 0; jj < n2; jj++) {
                    if (dat[ii][jj] != NULL) {
                        for (kk = 0; kk < n3; kk++) {
                            if (dat[ii][jj][kk] != NULL) {
                                for (ll = 0; ll < n4; ll++) {
                                    if (dat[ii][jj][kk][ll] != NULL) {
                                        for (mm = 0; mm < n5; mm++) {
                                            free(dat[ii][jj][kk][ll][mm]);
                                        }
                                        free(dat[ii][jj][kk][ll]);
                                    }
                                }
                                free(dat[ii][jj][kk]);
                            }
                        }
                        free(dat[ii][jj]);
                    }
                }
                free(dat[ii]);
            }

        }
    }
}

void free_xillTable(xillTable *tab) {
    if (tab != NULL) {

        int ii;
        if (tab->data_storage != NULL) {
            if (tab->num_param == 5) {
                free_xillTable_5dim(tab->data_storage[0],
                                    tab->num_param_vals[0],
                                    tab->num_param_vals[1],
                                    tab->num_param_vals[2],
                                    tab->num_param_vals[3],
                                    tab->num_param_vals[4]
                );
            } else {
                for (ii = 0; ii < tab->num_param_vals[0]; ii++) {
                    if (tab->data_storage[ii] != NULL) {
                        free_xillTable_5dim(tab->data_storage[ii],
                                            tab->num_param_vals[1],
                                            tab->num_param_vals[2],
                                            tab->num_param_vals[3],
                                            tab->num_param_vals[4],
                                            tab->num_param_vals[5]
                        );
                    }
                }
            }
            free(tab->data_storage);
        }

        if (tab->param_vals != NULL) {
            for (ii = 0; ii < tab->num_param; ii++) {
                free(tab->param_vals[ii]);
            }
        }

        free(tab->num_param_vals);

        free(tab);
    }
}


void free_cached_xillTable(void) {
    free_xillTable(cached_xill_tab);
    free_xillTable(cached_xill_tab_dens);
    free_xillTable(cached_xill_tab_nthcomp);
}
