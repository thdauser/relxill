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

#include "rellp.h"

lpTable* cached_lp_table = NULL;

/** routine for the broken power law emissivity **/
static void get_emis_bkn(double* emis, double* re,int nr,
		double index1, double index2, double rbr){

	// TODO: make it independent of rbr (maybe 1 at r=1rg?)
	double alpha;

	int ii;
	for (ii=0; ii<nr; ii++){
		alpha = index1;
		if (re[ii] > rbr){
			alpha = index2;
		}
		emis[ii] = pow(re[ii]/rbr,-alpha);
	}

	return;
}

/* function to check of the system parameters need to be re-calculated  */
static int redo_get_emis_lp(relParam* param, relParam* rel_param, relSysPar* sysPar){

	if (rel_param==NULL){
		return 1;
	}

	double cache_limit = 1e-8;
	if(rel_param!=NULL && param!=NULL){
		if (fabs(param->a - rel_param->a) > cache_limit) {
			return 1;
		}
		if (fabs(param->height - rel_param->height) > cache_limit) {
			return 1;
		}
		if (fabs(param->rin - rel_param->rin) > cache_limit) {
			return 1;
		}
		if (fabs(param->rout - rel_param->rout) > cache_limit) {
			return 1;
		}
		if (fabs(param->gamma - rel_param->gamma) > cache_limit) {
			return 1;
		}
	} else {
		return 1;
	}

	return 0;
}

static void get_emis_jet_point_source(relParam* param, double* emis, double* del_emit, double* del_inc,
		double* re, int n_r, lpTable* tab, int ind_a, double ifac_a, int* status){

	lpDat* dat[2];
	int ind_h[2];
	double ifac_h[2];

	// we set the
	int ii;
	for (ii=0; ii<2; ii++){
		dat[ii] = tab->dat[ind_a+ii];
		ind_h[ii] = binary_search_float(dat[ii]->h,tab->n_h,param->height);
		ifac_h[ii]   = (param->height-dat[ii]->h[ind_h[ii]])/
					   (dat[ii]->h[ind_h[ii]+1]-dat[ii]->h[ind_h[ii]]);
	}

	double jet_rad[tab->n_rad];
	double jet_emis[tab->n_rad];
	double jet_del[tab->n_rad];
	double jet_del_inc[tab->n_rad];

	// interpolate everything for the given a-h values on the original grid
	for (ii=0; ii<tab->n_rad; ii++){
	     // #1: intensity
	     jet_emis[ii]  =
	          (1.0-ifac_a)*(1.0-ifac_h[0])*dat[0]->intens[ind_h[0]][ii]
	          +(1.0-ifac_a)*(ifac_h[0])*dat[0]->intens[ind_h[0]+1][ii]
	          +(ifac_a)*(1.0-ifac_h[1])*dat[1]->intens[ind_h[1]][ii]
	          +(ifac_a)*(ifac_h[1])*dat[1]->intens[ind_h[1]+1][ii];

	     jet_del[ii]  =
	          (1.0-ifac_a)*(1.0-ifac_h[0])*dat[0]->del[ind_h[0]][ii]
	          +(1.0-ifac_a)*(ifac_h[0])*dat[0]->del[ind_h[0]+1][ii]
	          +(ifac_a)*(1.0-ifac_h[1])*dat[1]->del[ind_h[1]][ii]
	          +(ifac_a)*(ifac_h[1])*dat[1]->del[ind_h[1]+1][ii];

	     jet_del_inc[ii]  =
	          (1.0-ifac_a)*(1.0-ifac_h[0])*dat[0]->del_inc[ind_h[0]][ii]
	          +(1.0-ifac_a)*(ifac_h[0])*dat[0]->del_inc[ind_h[0]+1][ii]
	          +(ifac_a)*(1.0-ifac_h[1])*dat[1]->del_inc[ind_h[1]][ii]
	          +(ifac_a)*(ifac_h[1])*dat[1]->del_inc[ind_h[1]+1][ii];

	     // #2: r-grid
	     jet_rad[ii] = interp_lin_1d(ifac_a,dat[0]->rad[ii],dat[1]->rad[ii]);
	 }

	// and now rebin it to the given radial grid (TODO: check validity as we differ from the original code here)
	double inter_r;

	// get the extent of the disk (indices are defined such that tab->r[ind+1] <= r < tab->r[ind]
	int ind_rmin = binary_search(jet_rad,tab->n_rad,param->rin);
	assert(ind_rmin>0);
	int kk=ind_rmin;
	for (ii=n_r-1 ; ii>=0 ;ii--){
		while((re[ii] >= jet_rad[kk+1])){
			kk++;
			if (kk>=tab->n_rad-1) { //TODO: construct table such that we don't need this?
				if ( re[ii]-RELTABLE_MAX_R <= 1e-6){
					kk=tab->n_rad-2;
					break;
				} else {
					RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid",status);
					printf("   --> radius %.4e ABOVE the maximal possible radius of %.4e \n",
							re[ii], RELTABLE_MAX_R);
					CHECK_STATUS_VOID(*status);
				}
			}
		}


		// for larger angles logarithmic interpolation works slightly better
		if (jet_del[kk]/M_PI*180.0<=75.0){
			inter_r = (re[ii]-jet_rad[kk])/(jet_rad[kk+1]-jet_rad[kk]);
		} else {
			inter_r = (log(re[ii])-log(jet_rad[kk]))/
					(log(jet_rad[kk+1])-log(jet_rad[kk]));
		}

		//  log grid for the intensity (due to the function profile)
		emis[ii] = interp_log_1d(inter_r, jet_emis[kk], jet_emis[kk+1]);

		del_emit[ii] = interp_lin_1d(inter_r, jet_del[kk], jet_del[kk+1]);
		del_inc[ii] = interp_lin_1d(inter_r, jet_del_inc[kk], jet_del_inc[kk+1]);

	     // multiply by the additional factor gi^gamma (see Dauser et al., 2013)
		 // NEW: after Adam Ingram's comments we use the corect Lorentz invariant phase space density and get a +2
		 emis[ii] *= pow(gi_potential_lp(re[ii],param->a,param->height,param->beta,del_emit[ii]), param->gamma+2);


	     // take the beaming of the jet into account (see Dauser et al., 2013)
	     if (param->beta > 1e-6) {
	        emis[ii] *= pow(doppler_factor(del_emit[ii],param->beta),2);
	     }
	}

}

/** routine to calculate the emissivity in the lamp post geometry**/
void get_emis_jet(relParam* param, double* emis, double* del_emit, double* del_inc,
		double* re, int n_r, int* status){

	if (cached_lp_table==NULL){
		read_lp_table(LPTABLE_FILENAME, &cached_lp_table,status);
		CHECK_STATUS_VOID(*status);
	}

	// calculate the emissivity
	assert(emis!=NULL);

	int ind_a   = binary_search_float(cached_lp_table->a,cached_lp_table->n_a,param->a);
	double ifac_a   = (param->a-cached_lp_table->a[ind_a])/
				   (cached_lp_table->a[ind_a+1]-cached_lp_table->a[ind_a]);

	// TODO: make it possible for an extended jet
	get_emis_jet_point_source(param, emis, del_emit, del_inc, re, n_r,
			cached_lp_table, ind_a, ifac_a,status);

	return;
}

/** calculate the emissivity profile, depending on the EMIS_TYPE given **/
void calc_emis_profile(relParam* param, relParam* cached_param, relSysPar* sysPar, int* status){

	assert(sysPar!=NULL);
	double invalid_angle = -1.0;

	/**  *** Broken Power Law Emissivity ***  **/
	if (param->emis_type==EMIS_TYPE_BKN){

		get_emis_bkn(sysPar->emis, sysPar->re, sysPar->nr,
				param->emis1,param->emis2,param->rbr);
		// set the angles in this case to a default value
		int ii;
		for (ii=0; ii<sysPar->nr; ii++){
			sysPar->del_emit[ii] = invalid_angle;
			sysPar->del_inc[ii] = invalid_angle;
		}

	/**  *** Lamp Post Emissivity ***  **/
	} else if (param->emis_type==EMIS_TYPE_LP){

		if (redo_get_emis_lp(param, cached_param, sysPar)){
			get_emis_jet(param,sysPar->emis, sysPar->del_emit, sysPar->del_inc,
					sysPar->re, sysPar->nr, status);
		}

	} else {

		RELXILL_ERROR(" calculation of emissivity profile not possible \n",status);
		printf("   -> emis_type=%i not known \n",param->emis_type);
		return;
	}
	return;
}

void free_cached_lpTable(void){
	free_lpTable(cached_lp_table);
}
