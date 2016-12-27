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


/** global parameters, which can be used for several calls of the model */
relParam* cached_rel_param=NULL;
relTable* relline_table=NULL;
relSysPar* cached_relSysPar=NULL;
relSysPar* cached_tab_sysPar=NULL;
double cache_limit = 1e-8;

/* function to check of the system parameters need to be re-calculated  */
static int redo_get_system_parameters(relParam* param,relParam* cached_rel_param){

	if (cached_relSysPar==NULL){
		return 1;
	}

	if (cached_rel_param!=NULL){
		if (fabs(param->a - cached_rel_param->a) > cache_limit) {
			return 1;
		}

		if (fabs(param->incl - cached_rel_param->incl) > cache_limit) {
			return 1;
		}
	}

	return 0;
}

// interpolate the table in the A-MU0 plane (for one value of radius)
static void interpol_a_mu0(int ii, double ifac_a, double ifac_mu0, int ind_a,
		int ind_mu0, relSysPar* sysPar, relTable* relline_table) {
	sysPar->gmin[ii] = interp_lin_2d_float(ifac_a, ifac_mu0,
			relline_table->arr[ind_a][ind_mu0]->gmin[ii], relline_table->arr[ind_a+1][ind_mu0]->gmin[ii],
			relline_table->arr[ind_a][ind_mu0+1]->gmin[ii], relline_table->arr[ind_a+1][ind_mu0+1]->gmin[ii]);

	sysPar->gmax[ii] = interp_lin_2d_float(ifac_a, ifac_mu0,
			relline_table->arr[ind_a][ind_mu0]->gmax[ii], relline_table->arr[ind_a+1][ind_mu0]->gmax[ii],
			relline_table->arr[ind_a][ind_mu0+1]->gmax[ii], relline_table->arr[ind_a+1][ind_mu0+1]->gmax[ii]);

	for (int jj = 0; jj < relline_table->n_g; jj++) {
		sysPar->trff[ii][jj][0] = interp_lin_2d_float(ifac_a, ifac_mu0,
				relline_table->arr[ind_a][ind_mu0]->trff1[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0]->trff1[ii][jj],
				relline_table->arr[ind_a][ind_mu0+1]->trff1[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0+1]->trff1[ii][jj]);

		sysPar->trff[ii][jj][1] = interp_lin_2d_float(ifac_a, ifac_mu0,
				relline_table->arr[ind_a][ind_mu0]->trff2[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0]->trff2[ii][jj],
				relline_table->arr[ind_a][ind_mu0+1]->trff2[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0+1]->trff2[ii][jj]);

		sysPar->cosne[ii][jj][0] = interp_lin_2d_float(ifac_a, ifac_mu0,
				relline_table->arr[ind_a][ind_mu0]->cosne1[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0]->cosne1[ii][jj],
				relline_table->arr[ind_a][ind_mu0+1]->cosne1[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0+1]->cosne1[ii][jj]);

		sysPar->cosne[ii][jj][1] = interp_lin_2d_float(ifac_a, ifac_mu0,
				relline_table->arr[ind_a][ind_mu0]->cosne2[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0]->cosne2[ii][jj],
				relline_table->arr[ind_a][ind_mu0+1]->cosne2[ii][jj],
				relline_table->arr[ind_a+1][ind_mu0+1]->cosne2[ii][jj]);

	}
}

/*  get the fine radial grid */
static void get_fine_radial_grid(double rin, double rout, relSysPar* sysPar){

	double r1=1.0/sqrt(rout);
	double r2=1.0/sqrt(rin);
	for (int ii=0; ii<sysPar->nr; ii++){
		sysPar->re[ii] = ((double) (ii) )*(r2-r1)/sysPar->nr+r1;
		sysPar->re[ii] = 1.0/(sysPar->re[ii]*sysPar->re[ii]);
		assert(sysPar->re[ii]>1.0);
	}
	return;
}

/* function interpolating the rel table values for rin,rout,mu0,incl   */
static void interpol_relTable(relSysPar* sysPar,double a, double mu0, double rin, double rout,
		relTable* relline_table, int* status){

	// load tables
	if (relline_table==NULL){
		read_relline_table(RELTABLE_FILENAME,&relline_table,status);
		CHECK_STATUS_VOID(*status);

	}


	double rms = kerr_rms(a);

	// make sure the desired rmin is within bounds and order correctly
	assert(rout>rin);
	assert(rin>=rms);


	/**************************************/
	/** 1 **  Interpolate in A-MU0 plane **/
	/**************************************/

	// get a structure to store the values from the interpolation in the A-MU0-plane
	if (cached_tab_sysPar == NULL){
		cached_tab_sysPar = new_relSysPar(relline_table->n_r,relline_table->n_g,status);
		CHECK_STATUS_VOID(*status);
	}

	int ind_a   = binary_search_float(relline_table->a,relline_table->n_a,a);
	int ind_mu0 = binary_search_float(relline_table->mu0,relline_table->n_mu0,mu0);

	double ifac_a   = (a-relline_table->a[ind_a])/
				   (relline_table->a[ind_a+1]-relline_table->a[ind_a]);
	double ifac_mu0 = (mu0-relline_table->mu0[ind_mu0])/
				   (relline_table->mu0[ind_mu0+1]-relline_table->mu0[ind_mu0]);

	/** TODO: check if we need R_DIFF
	double r_diff = rms - ((1.0-ifac_a)*kerr_rms(relline_table->a[ind_a])
	       + (ifac_a)*kerr_rms(relline_table->a[ind_a+1]) ); **/

//	printf(" rms:%.3f, rdiff:%.3e -- %.3e\n",rms,r_diff,relline_table->arr[ind_a][ind_mu0]->r[relline_table->n_r-1]);
//	printf(" ind_a:%i ind_mu0:%i \n",ind_a,ind_mu0);

	/** get the radial grid (the radial grid only changes with A by the table definition) **/
	assert(fabs(relline_table->arr[ind_a][ind_mu0]->r[relline_table->n_r-1]
			    - relline_table->arr[ind_a][ind_mu0]->r[relline_table->n_r-1]) < 1e-6);
	for (int ii=0; ii < relline_table->n_r; ii++){
		cached_tab_sysPar->re[ii] = interp_lin_1d(ifac_a,
				relline_table->arr[ind_a][ind_mu0]->r[ii],relline_table->arr[ind_a+1][ind_mu0]->r[ii]);
	}

	// get the extent of the disk (indices are defined such that tab->r[ind] <= r < tab->r[ind+1]
	int ind_rmin = inv_binary_search(cached_tab_sysPar->re,relline_table->n_r,rin);
	int ind_rmax = inv_binary_search(cached_tab_sysPar->re,relline_table->n_r,rout);

	for (int ii=0; ii < relline_table->n_r; ii++){
		// TODO: SHOULD WE ONLY INTERPOLATE ONLY THE VALUES WE NEED??? //
		// only interpolate values where we need them
		if (ii>=ind_rmin || ii<=ind_rmax+1){
			interpol_a_mu0(ii, ifac_a, ifac_mu0, ind_a, ind_mu0, cached_tab_sysPar,relline_table);
		} else {  // set everything we won't need to 0 (just to be sure)
			cached_tab_sysPar->gmin[ii]=0.0;
			cached_tab_sysPar->gmax[ii]=0.0;
		    for (int jj=0; jj<relline_table->n_g;jj++){
			    for (int kk=0; kk<2;kk++){
			    	cached_tab_sysPar->trff[ii][jj][kk]=0.0;
			    	cached_tab_sysPar->cosne[ii][jj][kk]=0.0;
			    }
			 }
		}
	}

	/****************************/
	/** 2 **  Bin to Fine Grid **/
	/****************************/

	// only need to initialize and allocat memory if not already loaded
	if (sysPar == NULL){
		sysPar = new_relSysPar(N_FRAD,relline_table->n_g,status);
		CHECK_STATUS_VOID(*status);
	}
	get_fine_radial_grid(rin,rout,sysPar);

	// let's try to be as efficient as possible here (not that "r" DEcreases)
	assert(ind_rmin>0); // as defined inverse, re[ind_rmin+1] is the lowest value
	assert((cached_tab_sysPar->re[ind_rmin+1]<=rin));
	assert((cached_tab_sysPar->re[ind_rmin]>=rin));
	assert((cached_tab_sysPar->re[ind_rmax]<=rout));
	assert(ind_rmax <= ind_rmin);
	assert(rout<=1000.0);

	double ifac_r;
	int ind_tabr=ind_rmin;
	int ii;
	for (ii=sysPar->nr-1 ; ii>=0 ;ii--){
		while((sysPar->re[ii] >= cached_tab_sysPar->re[ind_tabr])){
			ind_tabr--;
			if (ind_tabr<0) { //TODO: construct table such that we don't need this?
				if ( sysPar->re[ii]-RELTABLE_MAX_R <= 1e-6){
					ind_tabr=0;
					break;
				} else {
					RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid",status);
					printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
							sysPar->re[ii], RELTABLE_MAX_R);
					CHECK_STATUS_VOID(*status);
				}
			}
		}

		ifac_r = (sysPar->re[ii]- cached_tab_sysPar->re[ind_tabr+1])
				/(cached_tab_sysPar->re[ind_tabr] - cached_tab_sysPar->re[ind_tabr+1]);
		assert(ifac_r>=0.0);

		// we only allow extrapolation (i.e. ifac_r < 0) for the last bin
		if (ifac_r >1.0 && ind_tabr>0){
			RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid",status);
			printf("   --> radius %.4e not found in [%.4e,%.4e]  \n",
					sysPar->re[ii],cached_tab_sysPar->re[ind_tabr+1],cached_tab_sysPar->re[ind_tabr]);
			CHECK_STATUS_VOID(*status);
		}

/**		printf("[%03i] rad: %.4e  <= %.4e  < %.4e [ifac_r=%.4e]\n",
				ii,cached_tab_sysPar->re[ind_tabr+1],sysPar->re[ii],
				cached_tab_sysPar->re[ind_tabr],ifac_r); **/

		for (int jj; jj<sysPar->ng; jj++){
			for (int kk; kk<2; kk++){
		        sysPar->trff[ii][jj][kk]=
		        	interp_lin_1d(ifac_r,cached_tab_sysPar->trff[ind_tabr+1][jj][kk]
										,cached_tab_sysPar->trff[ind_tabr][jj][kk]);

		        sysPar->cosne[ii][jj][kk]=
		        	interp_lin_1d(ifac_r,cached_tab_sysPar->cosne[ind_tabr+1][jj][kk]
										,cached_tab_sysPar->cosne[ind_tabr][jj][kk]);

			}
		}
 		sysPar->gmin[ii] =
 			interp_lin_1d(ifac_r,cached_tab_sysPar->gmin[ind_tabr+1],cached_tab_sysPar->gmin[ind_tabr]);
 		sysPar->gmax[ii] =
 			interp_lin_1d(ifac_r,cached_tab_sysPar->gmax[ind_tabr+1],cached_tab_sysPar->gmax[ind_tabr]);

	}

	return;
}

/* function to get the system parameters; decides if values need to be re-computed
 * interpolate values loaded from the table if necessary */
static relSysPar* get_system_parameters(relParam* param,relTable* relline_table, int* status){

	// only re-do the interpolation if rmin,rmax,a,mu0 changed
	// or if the cached parameters are NULL
	if (redo_get_system_parameters(param,cached_rel_param)){
		double mu0 = cos(param->incl);
		interpol_relTable(cached_relSysPar,param->a,mu0,param->rin,param->rout,relline_table,status);
		CHECK_STATUS_RET(*status,NULL);
	}

	return cached_relSysPar;
}


/* the relbase function calculating the basic relativistic line shape for a given parameter setup
 * input: ener(n_ener), param
 * output: photar(n_ener)     */
void relbase(const double* ener, const int n_ener, double* photar, relParam* param, int* status){

	// initialize parameter values
	relSysPar* sysPar = get_system_parameters(param,relline_table,status);
	CHECK_STATUS_VOID(*status);

	// get emissivity profile


	// get line shape


	// cache everything once we've calculated all the stuff
	// cached_params (!!!)
}

void free_cached_tables(void){
	free_relTable(relline_table);
	free_relSysPar(cached_relSysPar);
	free_relSysPar(cached_tab_sysPar);
}

relSysPar* new_relSysPar(int nr, int ng, int* status){
	relSysPar* sysPar = (relSysPar*) malloc( sizeof(relSysPar) );
	CHECK_MALLOC_RET_STATUS(sysPar,status,NULL);

	sysPar->ng = ng;
	sysPar->nr = nr;

	sysPar->re = (double*) malloc (nr*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->re,status,sysPar);
	sysPar->gmin = (double*) malloc (nr*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->gmin,status,sysPar);
	sysPar->gmax = (double*) malloc (nr*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->gmax,status,sysPar);

	sysPar->gstar = (double*) malloc (ng*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->gstar,status,sysPar);
	// we already set the values as they are fixed
	double H = 2e-3;
	for (int ii=0; ii<ng;ii++){
		sysPar->gstar[ii] = H + (1.0-2*H)/(ng-1)*( (float) (ii) );
	}

	sysPar->trff = (double***) malloc(nr*sizeof(double**));
	CHECK_MALLOC_RET_STATUS(sysPar->trff,status,sysPar);
	sysPar->cosne = (double***) malloc(nr*sizeof(double**));
	CHECK_MALLOC_RET_STATUS(sysPar->cosne,status,sysPar);
	for (int ii=0; ii < nr; ii++){
		sysPar->trff[ii] = (double**) malloc(ng*sizeof(double*));
		CHECK_MALLOC_RET_STATUS(sysPar->trff[ii],status,sysPar);
		sysPar->cosne[ii] = (double**) malloc(ng*sizeof(double*));
		CHECK_MALLOC_RET_STATUS(sysPar->cosne[ii],status,sysPar);
		for (int jj=0; jj < ng; jj++){
			sysPar->trff[ii][jj] = (double*) malloc(2*sizeof(double));
			CHECK_MALLOC_RET_STATUS(sysPar->trff[ii][jj],status,sysPar);
			sysPar->cosne[ii][jj] = (double*) malloc(2*sizeof(double));
			CHECK_MALLOC_RET_STATUS(sysPar->cosne[ii][jj],status,sysPar);
		}
	}

	return sysPar;
}

void free_relSysPar(relSysPar* sysPar){
	if (sysPar!=NULL){
		free(sysPar->re);
		free(sysPar->gmin);
		free(sysPar->gmax);
		free(sysPar->gstar);

		if(sysPar->trff != NULL){
			for (int ii=0; ii<sysPar->nr; ii++){
				if(sysPar->trff[ii] != NULL){
					for (int jj=0; jj<sysPar->ng; jj++){
						free(sysPar->trff[ii][jj]);
					}
					free(sysPar->trff[ii]);
				}
			}
			free(sysPar->trff);
		}

		if(sysPar->cosne != NULL){
			for (int ii=0; ii<sysPar->nr; ii++){
				if(sysPar->cosne[ii] != NULL){
					for (int jj=0; jj<sysPar->ng; jj++){
						free(sysPar->cosne[ii][jj]);
					}
					free(sysPar->cosne[ii]);
				}
			}
			free(sysPar->cosne);
		}

	}
}

