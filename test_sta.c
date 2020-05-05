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
#include "relbase.h"
#include "relutility.h"
#include "reltable.h"
#include "test_relxill.h"
#include "test_xilltab.c"

#define LIMIT_PREC 1e-6

/** test the currently implemented relline table
 ** [current version used: rel_table_v0.4e]   */
void test_relline_table(int* status){

	printf("\n ==> Testing RELLINE TABLE (%s) \n",RELTABLE_FILENAME);
	relTable* tab=NULL;

    // load the table
    read_relline_table(RELTABLE_FILENAME, &tab, status);
    CHECK_RELXILL_ERROR("loading the rel table failed", status);
    CHECK_STATUS_VOID(*status)

    // test certain values
    assert(tab != NULL);
    assert(tab->arr != NULL);

    double aref_val = 0.9829549;
    if (fabs(tab->a[RELTABLE_NA - 2] - aref_val) > LIMIT_PREC) {
        printf(" testing spin: expecting value of %f, but found %f\n",
               aref_val, tab->a[RELTABLE_NA - 2]);
        RELXILL_ERROR("values in rel table not correct", status);
    }

    double mu0ref_val = 0.09476821;
    if (fabs(tab->mu0[1] - mu0ref_val) > LIMIT_PREC) {
        printf(" testing mu0: expecting value of %f, but found %f\n", mu0ref_val, tab->mu0[1]);
        RELXILL_ERROR("values in rel table not correct", status);
    }

    const int n = 5;
    const float ref_val[5] = {(float) 985.76074, (float) 0.01127052,
                              (float) 0.01120779, (float) 0.01121218, (float) 0.03022655};
    const float val[5] = {
            tab->arr[0][0]->r[1],
            tab->arr[0][0]->trff1[0][0],
            tab->arr[0][0]->trff1[0][1],
            tab->arr[0][0]->trff1[1][0],
            tab->arr[0][1]->trff1[1][0]
    };
    int ii;
    for (ii = 0; ii < n; ii++) {
        if (fabsf(ref_val[ii] - val[ii]) > LIMIT_PREC) {
            printf(" testing rel table: expecting value of %f, but found %f\n", ref_val[ii], val[ii]);
            RELXILL_ERROR("values in rel table not correct", status);
            continue;
        }
    }
    printf("  ### all values of the relline table as expected\n");


    // free memory
    free_relTable(tab);

}


/** test the currently implemented relline table
 ** [current version used: rel_table_v0.4e]   */
void test_lp_table(int* status){

	printf("\n ==> Testing LP TABLE (%s) \n",LPTABLE_FILENAME);
	lpTable* tab=NULL;

    // load the table
    read_lp_table(LPTABLE_FILENAME, &tab, status);
    CHECK_RELXILL_ERROR("loading the rel table failed", status);

    // test certain values
    assert(tab != NULL);
    assert(tab->dat != NULL);
    assert(tab->a != NULL);

    double aref_val = 0.98374581;
    if (fabs(tab->a[LPTABLE_NA - 2] - aref_val) > LIMIT_PREC) {
        printf(" testing spin: expecting value of %f, but found %f\n",
               aref_val, tab->a[LPTABLE_NA - 2]);
        RELXILL_ERROR("values in rel table not correct", status);
    }

    double href_val = 1.8294963;
    if (fabs(tab->dat[1]->h[1] - href_val) > LIMIT_PREC) {
        printf(" testing hgrid: expecting value of %f, but found %f\n", href_val, tab->dat[1]->h[1]);
        RELXILL_ERROR("values in rel table not correct", status);
    }

    const int n = 3;
    const float ref_val[3] = {(float) 7.618106e-05, (float) 2.6826601, (float) -1.2509402};
    const float val[3] = {
            tab->dat[1]->intens[2][1],
            tab->dat[1]->del[2][1],
            tab->dat[1]->del_inc[2][1]
    };
    int ii;
    for (ii = 0; ii < n; ii++) {
        if (fabsf((ref_val[ii] - val[ii]) / ref_val[ii]) > LIMIT_PREC) {
            printf(" testing lp table: expecting value of %e, but found %e\n", ref_val[ii], val[ii]);
            RELXILL_ERROR("values in lp table not correct", status);
            continue;
        }
    }
    printf("  ### all values of the LP table as expected\n");

    // free memory
    free_lpTable(tab);

}


static void test_rebin_spectrum(int* status){

    int n0 = 6;
    int n = 4;
	double ener0[] = {1,3,4,5,6,7,9};
	double val0[] =  {1,1,1,2,1,1};

	double ener[] =  {0.5,2,4,5.5,7.5,8};
	double val[n];

	rebin_spectrum(ener,val,n,ener0,val0,n0);

	double val_ref[5];
	val_ref[0] = 0.5;
	val_ref[1] = 0.5+1;
	val_ref[2] = 1.0+1.0;
	val_ref[3] = 1.0+1.0+0.25;
	val_ref[4] = 0.25;

	printf("\n ==> Testing REBINNING Functions\n");

	int ii;
	for (ii=0; ii<n; ii++){
		if ( fabs(val[ii] - val_ref[ii] ) > LIMIT_PREC ){
			RELXILL_ERROR(" TEST-ERROR: testing of the function 'rebin_spectrum' failed",status);
            continue;
		}
	}

}

static void test_interp(int* status){

	printf("\n ==> Testing INTERPOLATION Routines \n");

	double ifac = 0.2;
	double rlo = 1.0;
	double rhi = 2.0;
	double val_1d = 1.2;
	if (fabs(interp_lin_1d(ifac,rlo,rhi)-val_1d) > LIMIT_PREC){
		RELXILL_ERROR(" TEST-ERROR: 1d interpolation does not produce the expected result",status);
	}


	double ifac1 = 0.4;
	double ifac2 = 0.8;
	double r11 = 1.0;
	double r12 = 2.0;
	double r21 = 2.0;
	double r22 = 4.0;
	double val_2d = 2.52;
	if (fabs(interp_lin_2d(ifac1,ifac2,r11,r12,r21,r22)-val_2d) > LIMIT_PREC){
		RELXILL_ERROR(" TEST-ERROR: 2d interpolation does not produce the expected result",status);
		printf(" VALUE=%e instead of %e \n",interp_lin_2d(ifac1,ifac2,r11,r12,r21,r22),val_2d);
	}

    float rf11 = (float) 1.0;
    float rf12 = (float) 2.0;
    float rf21 = (float) 2.0;
    float rf22 = (float) 4.0;
	if (fabs(interp_lin_2d_float(ifac1,ifac2,rf11,rf12,rf21,rf22)-val_2d) > LIMIT_PREC){
		RELXILL_ERROR(" TEST-ERROR: 2d interpolation (float) does not produce the expected result",status);
	}


}

static void do_std_test(int* status){

	test_relline_table(status);
    CHECK_STATUS_VOID(*status)
	printf("     ---> successful \n");

	test_lp_table(status);
    CHECK_STATUS_VOID(*status)
	printf("     ---> successful \n");


	test_interp(status);
    CHECK_STATUS_VOID(*status)
	printf("     ---> successful \n");

	test_rebin_spectrum(status);
    CHECK_STATUS_VOID(*status)
	printf("     ---> successful \n");

}


int main(int argc, char *argv[]){
	char *buf;
	int status = EXIT_SUCCESS;



	int do_all = 1;
	int do_relline = 0;
	int do_rellinelp = 0;
	int do_relxill = 0;
	int do_relxilllp = 0;
	int do_relxilllpion = 0;
	int do_relxilldens = 0;
	int do_relxillns = 0;
    int do_relxillco = 0;
	int do_relxilllpdens = 0;
	int do_relxillnthcomp= 0;
	int do_relxilllpnthcomp = 0;
	int do_relxilllpionnthcomp = 0;
	int do_relconv = 0;

	if (argc>=2){
		if (strcmp(argv[1],"version")==0){
			get_version_number(&buf,&status);
			printf("%s",buf);
			free(buf);
			return status;
		}

        if (strcmp(argv[1], "relline") == 0) {
            do_relline = 1;
            do_all = 0;
		} else if (strcmp(argv[1],"relconv")==0){
				do_relconv=1;
				do_all=0;
		} else if (strcmp(argv[1],"rellinelp")==0){
				do_rellinelp=1;
				do_all=0;
        } else if (strcmp(argv[1], "relxill") == 0) {
            do_relxill = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllp") == 0) {
            do_relxilllp = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilldens") == 0) {
            do_relxilldens = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxillNS") == 0) {
            do_relxillns = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxillCO") == 0) {
            do_relxillco = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpdens") == 0) {
            do_relxilllpdens = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxillCp") == 0) {
            do_relxillnthcomp = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpCp") == 0) {
            do_relxilllpnthcomp = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpion") == 0) {
            do_relxilllpion = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpionCp") == 0) {
            do_relxilllpionnthcomp = 1;
            do_all = 0;
        }


    }

	int n = 1;
	if (argc==3){
        n = (int) strtod(argv[2], NULL);
	}

	do{
		get_version_number(&buf,&status);
		printf("\n === Starting RELXILL Version %s === \n\n",buf);
		free(buf);

		if (do_all){
			do_std_test(&status);

			test_xilltables();
		}

		status=EXIT_SUCCESS;
		if (do_all || do_relline){
            status=EXIT_SUCCESS;
			std_eval_relline(&status,n);
			if (status==EXIT_SUCCESS) {
                printf("     ---> successful \n");
            }
		}

		if (do_all || do_rellinelp){
            status=EXIT_SUCCESS;
			std_eval_relline_lp(&status,1);
            if (status==EXIT_SUCCESS) {
                printf("     ---> successful \n");
            }
		}

        if (do_all || do_relxill){
            status=EXIT_SUCCESS;
            std_eval_relxill(&status,n);
            CHECK_STATUS_BREAK(status)
            if (status==EXIT_SUCCESS) {
                printf("     ---> successful \n");
            }
            bugtest_eval_relxill(&status);
            CHECK_STATUS_BREAK(status)
            if (status==EXIT_SUCCESS) {
                printf("     ---> successful \n");
            }
        }

        if (do_relconv){
			std_eval_relconv(&status,1);
            CHECK_STATUS_BREAK(status)
			std_eval_relconvlp(&status,1);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");

		}

		if (do_all) {
			std_eval_xillver(&status,1);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");

		}


		if (do_all || do_relxillns){
			std_eval_relxill_ns(&status,n);
			CHECK_STATUS_BREAK(status);
			printf("     ---> successful \n");
		}

        if (do_all || do_relxillco) {
            std_eval_relxill_co(&status, n);
            CHECK_STATUS_BREAK(status);
            printf("     ---> successful \n");
        }

        if (do_all || do_relxilllp) {
			std_eval_relxilllp(&status,n);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");
		}

		if (do_all || do_relxillnthcomp){
			std_eval_relxill_nthcomp(&status,n);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");
		}

		if (do_all || do_relxilllpnthcomp){
			std_eval_relxilllp_nthcomp(&status,n);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");
		}

		if (do_all || do_relxilldens){
			std_eval_relxilldens(&status,n);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");
		}

		if (do_all || do_relxilllpdens){
			std_eval_relxilllpdens(&status,n);
            CHECK_STATUS_BREAK(status)
            printf("     ---> successful \n");
        }

        if (do_all || do_relxilllpion) {
            std_eval_relxilllpion(&status, n);
            CHECK_STATUS_BREAK(status)
            printf("     ---> successful \n");

        }

        if (do_all || do_relxilllpionnthcomp) {
            std_eval_relxilllpion_nthcomp(&status, n);
            CHECK_STATUS_BREAK(status)
            printf("     ---> successful \n");

        }


        printf("\n ==> Cleaning up and freeing cached structures\n");
		free_cached_tables();
        free_cache();

	} while(0);

	if(status!=EXIT_SUCCESS){
		printf(" *** TESTING NOT SUCCESSFUL \n");
		// free tables
		printf( "\n *** Cleaning up and freeing cached structures\n");
		free_cached_tables();
		free_cache();
	}

  return status;
}
