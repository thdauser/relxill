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

    Copyright 2017 Thomas Dauser, Remeis Observatory & ECAP
*/
#include "relbase.h"
#include "relmodels.h"
#include "relutility.h"
#include "reltable.h"


#define LIMIT_PREC 1e-6

static void set_std_param_xillver(double* inp_par, int* status){
	inp_par[0] = 2.0;
	inp_par[1] = 1.0;
	inp_par[2] = 1.0;
	inp_par[3] = 100.0;
	inp_par[4] = 45.0;
	inp_par[5] = 0.;
}


static void set_std_param_relline(double* inp_par, int* status){
	inp_par[0] = 1.0;
	inp_par[1] = 3.0;
	inp_par[2] = 3.0;
	inp_par[3] = 15.0;
	inp_par[4] = 0.998;
	inp_par[5] = 30.0;
	inp_par[6] = -1.;
	inp_par[7] = 400.;
	inp_par[8] = 0.0;
}


static void set_std_param_relline_lp(double* inp_par, int* status){
	inp_par[0] = 1.0;
	inp_par[1] = 3.0;
	inp_par[2] = 0.998;
	inp_par[3] = 30.0;
	inp_par[4] = -1.;
	inp_par[5] = 400.;
	inp_par[6] = 0.0;
	inp_par[7] = 2.0;
}

/** standard evaluation of the relline model **/
static void std_eval_relline(int* status, int n){

	printf("\n *** Evaluating RELLINE MODEL \n");
	/* set the parameters */
	int n_param = NUM_PARAM_RELLINE;
	double inp_par[NUM_PARAM_RELLINE];
	set_std_param_relline(inp_par, status);
	CHECK_STATUS_VOID(*status);

	/* create an energy grid */
	int n_ener = 300;
	double ener[n_ener+1];
	get_log_grid(ener,n_ener+1,0.1,1.5);

	/* call the relline model */
	double photar[n_ener];
	int ii;
	for (ii=0; ii<n; ii++){
		relline(ener,n_ener,photar,inp_par,n_param,status);
	}
}

/** standard evaluation of the relline model **/
static void std_eval_relline_lp(int* status, int n){

	printf("\n *** Evaluating RELLINE LP MODEL \n");
	/* set the parameters */
	int n_param = NUM_PARAM_RELLINELP;
	double inp_par_lp[n_param];
	set_std_param_relline_lp(inp_par_lp, status);
	CHECK_STATUS_VOID(*status);

	/* create an energy grid */
	int n_ener = 600;
	double ener[n_ener+1];
	get_log_grid(ener,n_ener+1,0.1,1.5);

	/* call the relline model */
	double photar[n_ener];
	int ii;
	for (ii=0; ii<n; ii++){
		rellinelp(ener,n_ener,photar,inp_par_lp,n_param,status);
		CHECK_STATUS_VOID(*status);
	}


}




/** standard evaluation of the relline model **/
static void std_eval_xillver(int* status, int n){

	printf("\n *** Evaluating XILLVER MODEL \n");
	/* set the parameters */
	int n_param = NUM_PARAM_XILLVER;
	double inp_par[n_param];
	set_std_param_xillver(inp_par, status);
	CHECK_STATUS_VOID(*status);

	/* create an energy grid */
	int n_ener = 1500;
	double ener[n_ener+1];
	get_log_grid(ener,n_ener+1,0.05,100.0);

	/* call the relline model */
	double photar[n_ener];
	int ii;
	for (ii=0; ii<n; ii++){
		xillver(ener,n_ener,photar,inp_par,n_param,status);
		CHECK_STATUS_VOID(*status);
	}


}


/** test the currently implemented relline table
 ** [current version used: rel_table_v0.4e]   */
void test_relline_table(int* status){

	printf("\n *** Testing RELLINE TABLE (%s) \n",RELTABLE_FILENAME);
	relTable* tab=NULL;

	do{
		// load the table
		read_relline_table(RELTABLE_FILENAME,&tab,status);
		CHECK_RELXILL_ERROR("loading the rel table failed",status);

		// test certain values
		assert(tab!=NULL);
		assert(tab->arr!=NULL);

		double aref_val = 0.98605162;
		if ( fabs(tab->a[RELTABLE_NA-2] - aref_val) > LIMIT_PREC ){
			printf(" testing spin: expecting value of %f, but found %f\n",
					aref_val,tab->a[RELTABLE_NA-2]);
			RELXILL_ERROR("values in rel table not correct",status);
			break;
		}

		double mu0ref_val = 9.7182781e-02;
		if ( fabs(tab->mu0[1] - mu0ref_val) > LIMIT_PREC ){
			printf(" testing mu0: expecting value of %f, but found %f\n",mu0ref_val,tab->mu0[1]);
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
		printf("  ### all values of the relline table as expected\n");

	} while(0);

	// free memory
	free_relTable(tab);

	return;
}


/** test the currently implemented relline table
 ** [current version used: rel_table_v0.4e]   */
void test_lp_table(int* status){

	printf("\n *** Testing LP TABLE (%s) \n",LPTABLE_FILENAME);
	lpTable* tab=NULL;

	do{
		// load the table
		read_lp_table(LPTABLE_FILENAME,&tab,status);
		CHECK_RELXILL_ERROR("loading the rel table failed",status);

		// test certain values
		assert(tab!=NULL);
		assert(tab->dat!=NULL);
		assert(tab->a!=NULL);

		double aref_val = 0.98374581;
		if ( fabs(tab->a[LPTABLE_NA-2] - aref_val) > LIMIT_PREC ){
			printf(" testing spin: expecting value of %f, but found %f\n",
					aref_val,tab->a[LPTABLE_NA-2]);
			RELXILL_ERROR("values in rel table not correct",status);
			break;
		}

		double href_val = 1.8294963;
		if ( fabs(tab->dat[1]->h[1] - href_val) > LIMIT_PREC ){
			printf(" testing hgrid: expecting value of %f, but found %f\n",href_val,tab->dat[1]->h[1]);
			RELXILL_ERROR("values in rel table not correct",status);
			break;
		}

		const int n = 3;
		const float ref_val[3] = {7.5840981e-05, 2.6826601, -1.2509402 };
		const float val[3] = {
				tab->dat[1]->intens[2][1],
				tab->dat[1]->del[2][1],
				tab->dat[1]->del_inc[2][1]
			};
		for (int ii=0; ii<n; ii++){
			if ( fabs( (ref_val[ii] - val[ii]) / ref_val[ii] ) > LIMIT_PREC ){
				printf(" testing lp table: expecting value of %f, but found %f\n",ref_val[ii],val[ii]);
				RELXILL_ERROR("values in lp table not correct",status);
				break;
			}
		}
		printf("  ### all values of the LP table as expected\n");

	} while(0);

	// free memory
	free_lpTable(tab);

	return;
}


static void test_interp(int* status){

	printf("\n *** Testing INTERPOLATION Routines \n");

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

	float rf11 = 1.0;
	float rf12 = 2.0;
	float rf21 = 2.0;
	float rf22 = 4.0;
	if (fabs(interp_lin_2d_float(ifac1,ifac2,rf11,rf12,rf21,rf22)-val_2d) > LIMIT_PREC){
		RELXILL_ERROR(" TEST-ERROR: 2d interpolation (float) does not produce the expected result",status);
	}


}

int main(void){
	char *buf;
	int status = EXIT_SUCCESS;

	do{
		get_version_number(&buf,&status);
		printf("\n *** Starting RELXILL Version %s *** \n\n",buf);
		free(buf);

//		test_relline_table(&status);
//		CHECK_STATUS_BREAK(status);
//		printf("     ---> successful \n");

/*		test_lp_table(&status);
		CHECK_STATUS_BREAK(status);
		printf("     ---> successful \n"); */

/*
		test_interp(&status);
		CHECK_STATUS_BREAK(status);
		printf("     ---> successful \n"); */


/*		std_eval_relline(&status,1);
		CHECK_STATUS_BREAK(status);
		printf("     ---> successful \n"); */

/*		std_eval_relline_lp(&status,1);
		CHECK_STATUS_BREAK(status);
		printf("     ---> successful \n"); */

		std_eval_xillver(&status,1);
		CHECK_STATUS_BREAK(status);
		printf("     ---> successful \n");


		printf( "\n *** Cleaning up and freeing cached structures\n");
		free_cached_tables();

	} while(0);

	if(status!=EXIT_SUCCESS){
		printf(" *** TESTING NOT SUCCESSFUL \n");
		// free tables
		printf( "\n *** Cleaning up and freeing cached structures\n");
		free_cached_tables();
	}

}
