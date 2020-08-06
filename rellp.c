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


/*
 * ignores param->height and param->beta
 */
static void get_emis_jet_point_source(relParam* param, double height, double beta, double* emis, double* del_emit, double* del_inc,
		double* re, int n_r, lpTable* tab, int ind_a, double ifac_a, int* status){

	lpDat* dat[2];
	int ind_h[2];
	double ifac_h[2];

	// we set the
	int ii;
	for (ii=0; ii<2; ii++){
		dat[ii] = tab->dat[ind_a+ii];
		ind_h[ii] = binary_search_float(dat[ii]->h,tab->n_h, (float) height);
		ifac_h[ii]   = (height-dat[ii]->h[ind_h[ii]])/
					   (dat[ii]->h[ind_h[ii]+1]-dat[ii]->h[ind_h[ii]]);

		// make sure the incident angle is defined as positive value (otherwise the interpolation
		// will create problems / jumps )
		int jj; int kk;
		for (jj=0; jj<tab->n_h; jj++){
			for (kk=0; kk<tab->n_rad; kk++){
				dat[ii]->del_inc[jj][kk] = fabs(dat[ii]->del_inc[jj][kk]);
				dat[ii]->del[jj][kk] = fabs(dat[ii]->del[jj][kk]);
			}
		}

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
//	int ind_rmin = binary_search(jet_rad,tab->n_rad,param->rin);
	int ind_rmin = binary_search(jet_rad,tab->n_rad,re[n_r-1]);
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

	    /** multiply by the additional factor gi^gamma (see Dauser et al., 2013) **/
		 emis[ii] *= pow(gi_potential_lp(re[ii],param->a,height,beta,del_emit[ii]), param->gamma);

	     // take the beaming of the jet into account (see Dauser et al., 2013)
	     if (param->beta > 1e-6) {
	        emis[ii] *= pow(doppler_factor(del_emit[ii],beta),2);
	     }
	}

}

static int modelPointsource(relParam* param){
  double htopPrecLimit = 1e-3;
  if (fabs(param->htop ) <=1.0 || ( param->htop+htopPrecLimit <= param->height)  ){
    return 1;
  } else {
    return 0;
  }
}


/*
 * Calculate the velocity (in units of c) for a given height, hbase and a given
 * velocity (beta) at 100Rg. Specfici definitions are:
 *  - beta is defined as the velocity at 100r_g above the jet base at hbase
 *  - we assume constant acceleration such that currentBeta=bet100
 */
double jetSpeedConstantAccel(double beta100, double height, double hbase){


  const double HEIGHT_REF = 100.0;  // in units of Rg

  double rel_gamma = 1.0/sqrt(1-pow(beta100,2));
  double acceleration = 1.0/HEIGHT_REF* (sqrt(1 + pow(rel_gamma*beta100, 2))  - 1  );

  double x = height-hbase;
  double t = sqrt(pow(x,2) + 2*x/acceleration);

  // speed for constant accel in SRT
  double beta = acceleration*t / sqrt(1 + (pow(acceleration*t, 2)));

  return beta;
}


/*
 *  The extended emission is calculated for a certain set of parameters / definitions
 *  - if htop <= heigh=hbase we assume it's a point-like jet
 *  - the meaning of beta for the extended jet is the velocity at 100Rg, in case the
 *    profile is of interest, it will be output in the debug mode
 */
void get_emis_jet_extended(relParam* param, double* emis, double* del_emit, double* del_inc,
                                      double* re, int n_r, lpTable* tab, int ind_a, double ifac_a, int* status) {


  // check and set the parameters as defined for the extended jet
  assert(param->height < param->htop);
  double beta100Rg = param->beta;
  double hbase = param->height;

  int nh = 100;
  double heightArray[nh+1];
  get_log_grid(heightArray, nh, hbase, param->htop);

  double singleEmis[n_r];

  for (int jj; jj<n_r; jj++){
    emis[jj] = 0.0;
  }

  double height[nh];
  double beta[nh];

  for (int ii = 0; ii < nh; ii++){

    height[ii] = 0.5*(heightArray[ii]+heightArray[ii+1]);
    beta[ii] = jetSpeedConstantAccel(beta100Rg, height[ii], hbase);

    get_emis_jet_point_source(param, height[ii], beta[ii], singleEmis, del_emit, del_inc, re, n_r,
                              cached_lp_table, ind_a, ifac_a, status);

    // assuming an constant luminosity in the frame of the jet
    double heightSegmentIntegrationFactor = (heightArray[ii+1] - heightArray[ii]) / (param->htop - hbase);

    for (int jj; jj<n_r; jj++){
      emis[jj] = singleEmis[jj]*heightSegmentIntegrationFactor;
    }


  }


  if (is_debug_run() ){
    save_radial_profile("test_rellxill_heightVelocityProfile.txt", height, beta, nh);
  }

}


/** routine to calculate the emissivity in the lamp post geometry**/
void get_emis_jet(relParam* param, double* emis, double* del_emit, double* del_inc,
		double* re, int n_r, int* status){

	CHECK_STATUS_VOID(*status);

	if (cached_lp_table==NULL){
		read_lp_table(LPTABLE_FILENAME, &cached_lp_table,status);
		CHECK_STATUS_VOID(*status);
	}

	// calculate the emissivity
	assert(emis!=NULL);

	int ind_a   = binary_search_float(cached_lp_table->a,cached_lp_table->n_a,(float) param->a);
	double ifac_a   = (param->a-cached_lp_table->a[ind_a])/
				   (cached_lp_table->a[ind_a+1]-cached_lp_table->a[ind_a]);


	if (modelPointsource(param)) {
      get_emis_jet_point_source(param, param->height, param->beta, emis, del_emit, del_inc, re, n_r,
                                cached_lp_table, ind_a, ifac_a, status);
    } else {
	  get_emis_jet_extended(param, emis, del_emit, del_inc, re, n_r,
                            cached_lp_table, ind_a, ifac_a, status);
	}

}


// calculate the angles of emission from the primary source to git Rin and Rout
void get_ad_del_lim(relParam* param, relSysPar* sysPar, int* status) {
	int nr_lim = 2;
	double del_emit[2];
	double del_dummy[2];
	double emis_dummy[2];
	double rad[2] = { AD_ROUT_MAX, kerr_rms(param->a) };
	// get the primary source emission angle for the simulated inner and out edge of the disk
	get_emis_jet(param, emis_dummy, del_emit, del_dummy, rad, nr_lim, status);
	sysPar->del_ad_rmax = del_emit[0];
	sysPar->del_ad_risco = del_emit[1];

	if (*status!=EXIT_SUCCESS){
		printf(" *** failed calculating the primary source photon emission angles for Rin=%.2e and Rmax=%.2e",rad[1], rad[0]);
	}
}


emisProfile* calc_emis_profile(double* re, int nr, relParam* param, int* status){

	CHECK_STATUS_RET(*status, NULL);

	emisProfile* emis = new_emisProfile(re, nr, status);

	double invalid_angle = -1.0;

	/**  *** Broken Power Law Emissivity ***  **/
	if (param->emis_type==EMIS_TYPE_BKN){

		get_emis_bkn(emis->emis, emis->re, emis->nr,
				param->emis1,param->emis2,param->rbr);
		// set the angles in this case to a default value
		int ii;
		for (ii=0; ii<nr; ii++){
			emis->del_emit[ii] = invalid_angle;
			emis->del_inc[ii] = invalid_angle;
		}

	/**  *** Lamp Post Emissivity ***  **/
	} else if (param->emis_type==EMIS_TYPE_LP){

      // TODO: it does not work for an extended jet to calculate del_emis and del_inc
		get_emis_jet(param,emis->emis, emis->del_emit, emis->del_inc,
				re, nr, status);
	} else {

		RELXILL_ERROR(" calculation of emissivity profile not possible \n",status);
		printf("   -> emis_type=%i not known \n",param->emis_type);
		return NULL;
	}


	if (*status!=EXIT_SUCCESS){
		RELXILL_ERROR("calculating the emissivity profile failed due to previous error", status);
	}

	return emis;
}

void free_cached_lpTable(void){
	free_lpTable(cached_lp_table);
}
