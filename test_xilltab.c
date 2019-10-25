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
#include "relbase.h"
#include "relutility.h"
#include "test_relxill.h"


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

static void test_init_xilltable(char *fname, xillParam *param, int *status) {

    printf("\n *** TEST: initializing xilltable %s  ...   \n\n", fname);

    xillTable *tab = NULL;
    init_xillver_table(fname, &tab, param, status);

    assert(tab->n_ener > 0);

    if (tab->num_param == 5) {
        assert(tab->dat[0][0][0][0][0][0] == NULL);
    }

    if (*status != EXIT_SUCCESS) {
        printf("  FAILED");
    } else {
        printf("  SUCCESSFUL");
    }

}

static void test_all_xilltables(int *status) {

    xillParam *param = get_std_param_xillver(status);
    test_init_xilltable(XILLTABLE_FILENAME, param, status);

}

int main(int argc, char *argv[]) {
    char *buf;
    int status = EXIT_SUCCESS;

    get_version_number(&buf, &status);
    printf("\n === Testing XILLVER with RELXILL Version %s === \n\n", buf);
    free(buf);

    test_new_xilltable(&status);

    test_all_xilltables(&status);


    if (status != EXIT_SUCCESS) {
        printf(" *** TESTING NOT SUCCESSFUL \n");
        printf("\n *** Cleaning up and freeing cached structures\n");
    } else {
        printf(" *** TESTING SUCCESSFUL \n");
    }

}
