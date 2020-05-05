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
#include "test_relxill.h"
#include "relbase.h"


static void test_new_xilltable(int *status) {

    int dim[] = {5, 6};
    int ndim = 2;

    int ii;
    for (ii = 0; ii < ndim; ii++) {
        xillTable *tab = new_xillTable(dim[ii], status);

        free_xillTable(tab);
    }

    if (*status == EXIT_SUCCESS) {
        printf("\n *** TEST: setting up new XILLVER Table structure successful \n");
    } else {
        printf("\n *** TEST ERROR *** setting up new XILLVER Table structure NOT successful \n");
    }

}

static void print_status_outcome(const int *status) {
    if (*status != EXIT_SUCCESS) {
        printf("  FAILED \n\n");
    } else {
        printf("  SUCCESSFUL \n\n");
    }
}

static void test_tab_num_param(const xillParam *param, int *status, const xillTable *tab) {
    printf("     - number of parameters: %i ", tab->num_param);
    int dim_expect = 5;
    if (is_6dim_table(param->model_type)) {
        dim_expect = 6;
    }
    printf(" [expected: %i] ", dim_expect);
    if (tab->num_param != dim_expect) {
        *status = EXIT_FAILURE;
    }
}

static void test_init_xilltable(char *fname, xillParam *param, int *status) {


    xillTable *tab = NULL;
    init_xillver_table(fname, &tab, param, status);

    printf("\n *** TEST: initializing xilltable %s  ...  ", fname);

    if (*status==EXIT_SUCCESS) {
      assert(tab->n_ener > 0);

      int ii;
      for (ii = 0; ii < tab->num_param; ii++) {
        assert(tab->num_param_vals[ii] > 0);
      }

      assert(tab->num_param > 0);
      test_tab_num_param(param, status, tab);
    }

    print_status_outcome(status);

}

static void test_spec_norm(xill_spec *spec, int *status) {

    printf("\n *** TEST: check Spectrum normalization ...  ");

    int ii;
    int jj;

    double sum = 0.0;

    for (jj = 0; jj < spec->n_incl; jj++) {
        for (ii = 0; ii < spec->n_ener; ii++) {
            sum += spec->flu[jj][ii];
        }
    }

    printf("SUM = %.4e \n\n", sum);

}

static void test_get_spec(int *status, xillParam *param) {

    printf("\n *** TEST: loading Spectrum from Storage (model type = %i) ...   ",
           param->model_type);

    xill_spec *spec = get_xillver_spectra(param, status);
    if (*status != EXIT_SUCCESS) {
        printf("\n *** ERROR trying to load the spectrum \n");
    }

    assert(spec->n_ener > 0);
    assert(spec->n_incl > 0);

    print_status_outcome(status);


    test_spec_norm(spec, status);


    free_xill_spec(spec);
}

static void test_all_spec(int *status) {

    xillParam *param;

    // param = get_std_param_xillver(status);
    //  test_get_spec(status, param);

    param = get_std_param_xillver_co(status);
    test_get_spec(status, param);


}

static void test_all_xilltables(int *status) {

    xillParam *param;
    putenv("DEBUG_RELXILL=1");

    /** -1- xillver (5dim table) **/
    param = get_std_param_xillver(status);
    test_init_xilltable(XILLTABLE_FILENAME, param, status);

    /** -2- xillver CO  (6dim table) **/
    param = get_std_param_xillver_co(status);
    test_init_xilltable(XILLTABLE_CO_FILENAME, param, status);


    putenv("DEBUG_RELXILL=0");

}

void test_xilltables(void ) {
    char *buf;
    int status = EXIT_SUCCESS;

    get_version_number(&buf, &status);
    printf("\n === Testing XILLVER with RELXILL Version %s === \n\n", buf);
    free(buf);

    test_new_xilltable(&status);

    test_all_xilltables(&status);

    test_all_spec(&status);

    if (status != EXIT_SUCCESS) {
        printf(" *** TESTING XILLVER TABLES NOT SUCCESSFUL \n");
        printf("\n *** Cleaning up and freeing cached structures\n");
    } else {
        printf(" *** TESTING XILLVER TABLES SUCCESSFUL \n");
    }

}
