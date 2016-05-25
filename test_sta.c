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

    Copyright 2016 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "relbase.h"
#include "relmodels.h"
#include "relutility.h"
#include "reltable.h"

#define LIMIT_PREC 1e-8

static void set_std_param_relline(double* inp_par, int* status){
	inp_par[0] = 6.4;
	inp_par[1] = 3.0;
	inp_par[2] = 3.0;
	inp_par[3] = 5.0;
	inp_par[4] = 0.9;
	inp_par[5] = 30.;
	inp_par[6] = -1.;
	inp_par[7] = 400.;
	inp_par[8] = 0.0;
}

/** standard evaluation of the relline model **/
static void std_eval_relline(int* status){

	printf("\n *** Testing RELLINE MODEL \n");
	/* set the parameters */
	int n_param = NUM_PARAM_RELLINE;
	double inp_par[NUM_PARAM_RELLINE];
	set_std_param_relline(inp_par, status);
	CHECK_STATUS_VOID(*status);

	/* create an energy grid */
	int n_ener = 1000;
	double ener[n_ener];
	get_log_grid(ener,n_ener,0.1,10);

	/* call the relline model */
	double photar[n_ener];
	relline(ener,n_ener,photar,inp_par,n_param,status);

}

/** test the currently implemented relline table*/
void test_relline_table(int* status){

	printf("\n *** Testing RELLINE TABLE (%s) \n",RELTABLE_FILENAME);
	relTable* tab=NULL;

	do{
		// load the table
		read_relline_table(RELTABLE_FILENAME,&tab,status);
		CHECK_RELXILL_ERROR("loading the rel table failed",status);

		// test certain values
		assert(tab->arr!=NULL);

		double aref_val = 0.98605162;
		if ( abs(tab->a[RELTABLE_NA-2] - aref_val) > LIMIT_PREC ){
			printf(" testing spin: expecting value of %f, but found %f\n",aref_val,tab->a[RELTABLE_NA-1]);
			RELXILL_ERROR("values in rel table not correct",status);
			break;
		}

		double mu0ref_val = 0.98480862;
		if ( abs(tab->mu0[1] - mu0ref_val) > LIMIT_PREC ){
			printf(" testing mu0: expecting value of %f, but found %f\n",mu0ref_val,tab->mu0[RELTABLE_NA-1]);
			RELXILL_ERROR("values in rel table not correct",status);
			break;
		}

		const int n = 5;
		const float ref_val[5] = {985.76074, 0.0055350405,0.0055231834, 0.0054442692,0.030979944 };
		const float val[5] = {
				tab->arr[0][0]->r[1],
				tab->arr[0][0]->trff1[0][0],
				tab->arr[0][0]->trff1[0][1],
				tab->arr[0][0]->trff1[1][0],
				tab->arr[0][1]->trff1[1][0]
			};
		for (int ii=0; ii<n; ii++){
			if ( fabs( ref_val[ii] - val[ii]) > LIMIT_PREC ){
				printf(" testing rel table: expecting value of %f, but found %f\n",ref_val[ii],val[ii]);
				RELXILL_ERROR("values in rel table not correct",status);
				break;
			}
		}

	} while(0);

	// free memory
	free_relTable(tab);

	return;
}

int main(void){
	char *buf;
	int status = EXIT_SUCCESS;

	do{
		get_version_number(&buf,&status);
		printf("\n *** Starting RELXILL Version %s *** \n\n",buf);
		free(buf);

		test_relline_table(&status);
		CHECK_STATUS_BREAK(status);
		printf("     ---> successful \n");


		std_eval_relline(&status);
		CHECK_STATUS_BREAK(status);
		printf("     ---> successful \n");

	} while(0);

	if(status!=EXIT_SUCCESS){
		printf(" *** TESTING NOT SUCCESSFUL \n");
	}

}
