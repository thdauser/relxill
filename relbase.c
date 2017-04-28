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



/** global parameters, which can be used for several calls of the model */
relParam* cached_rel_param=NULL;
relTable* relline_table=NULL;
relSysPar* cached_relSysPar=NULL;
relSysPar* cached_tab_sysPar=NULL;
rel_spec* cached_rel_spec=NULL;

int save_1eV_pos = 0;
double cached_int_romb_rad = -1.0;
const double cache_limit = 1e-8;

const double ener_xill_norm_lo = 0.1;
const double ener_xill_norm_hi = 1000;
const int n_ener_xill = 3000;
double* ener_xill = NULL;

const int n_ener_std = N_ENER_CONV;
double* ener_std = NULL;


// precision to calculate gstar from [H:1-H] instead of [0:1]
const double GFAC_H = 2e-3;


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
		if (fabs(param->rin - cached_rel_param->rin) > cache_limit) {
			return 1;
		}
		if (fabs(param->rout - cached_rel_param->rout) > cache_limit) {
			return 1;
		}
		if (fabs(param->rbr - cached_rel_param->rbr) > cache_limit) {
			return 1;
		}
		if (fabs(param->height - cached_rel_param->height) > cache_limit) {
			return 1;
		}
		if (fabs(param->emis1 - cached_rel_param->emis1) > cache_limit) {
			return 1;
		}
		if (fabs(param->emis2 - cached_rel_param->emis2) > cache_limit) {
			return 1;
		}
	} else {
		return 1;
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

	int jj;
	for (jj = 0; jj < relline_table->n_g; jj++) {
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

	double r1=1.0/rout;
	double r2=1.0/rin;
	int ii;
	for (ii=0; ii<sysPar->nr; ii++){
		sysPar->re[ii] = ((double) (ii) )*(r2-r1)/(sysPar->nr-1)+r1;
		sysPar->re[ii] = 1.0/(sysPar->re[ii]);
		assert(sysPar->re[ii]>1.0);
	}
	return;
}

/* function interpolating the rel table values for rin,rout,mu0,incl   */
static void interpol_relTable(relSysPar** sysPar_inp,double a, double mu0, double rin, double rout,
		 int* status){

	// load tables
	if (relline_table==NULL){
		print_version_number(status);
		CHECK_STATUS_VOID(*status);
		read_relline_table(RELTABLE_FILENAME,&relline_table,status);
		CHECK_STATUS_VOID(*status);
	}
	relTable* tab = relline_table;
	assert(tab!=NULL);


	double rms = kerr_rms(a);

	// make sure the desired rmin is within bounds and order correctly
	assert(rout>rin);
	assert(rin>=rms);


	/**************************************/
	/** 1 **  Interpolate in A-MU0 plane **/
	/**************************************/

	// get a structure to store the values from the interpolation in the A-MU0-plane
	if (cached_tab_sysPar == NULL){
		cached_tab_sysPar = new_relSysPar(tab->n_r,tab->n_g,status);
		CHECK_STATUS_VOID(*status);
	}

	int ind_a   = binary_search_float(tab->a,tab->n_a,a);
	int ind_mu0 = binary_search_float(tab->mu0,tab->n_mu0,mu0);

	double ifac_a   = (a-tab->a[ind_a])/
				   (tab->a[ind_a+1]-tab->a[ind_a]);
	double ifac_mu0 = (mu0-tab->mu0[ind_mu0])/
				   (tab->mu0[ind_mu0+1]-tab->mu0[ind_mu0]);

	/** get the radial grid (the radial grid only changes with A by the table definition) **/
	assert(fabs(tab->arr[ind_a][ind_mu0]->r[tab->n_r-1]
			    - tab->arr[ind_a][ind_mu0]->r[tab->n_r-1]) < 1e-6);
	int ii;
	for (ii=0; ii < tab->n_r; ii++){
		cached_tab_sysPar->re[ii] = interp_lin_1d(ifac_a,
				tab->arr[ind_a][ind_mu0]->r[ii],tab->arr[ind_a+1][ind_mu0]->r[ii]);
	}
	// we have problems for the intermost radius due to linear interpolation (-> set to RISCO)
	if ( (cached_tab_sysPar->re[tab->n_r-1] > kerr_rms(a)) &&
			 ( (cached_tab_sysPar->re[tab->n_r-1] - kerr_rms(a)) / cached_tab_sysPar->re[tab->n_r-1]  < 1e-3)) {
		//		printf(" re-setting RIN from %.3f to %.3f (dr = %.2e)\n",cached_tab_sysPar->re[tab->n_r-1],kerr_rms(a),
		//		(cached_tab_sysPar->re[tab->n_r-1]-kerr_rms(a))/cached_tab_sysPar->re[tab->n_r-1]);
		cached_tab_sysPar->re[tab->n_r-1] = kerr_rms(a);
	}

	// get the extent of the disk (indices are defined such that tab->r[ind+1] <= r < tab->r[ind]
	int ind_rmin = inv_binary_search(cached_tab_sysPar->re,tab->n_r,rin);
	int ind_rmax = inv_binary_search(cached_tab_sysPar->re,tab->n_r,rout);

	int jj; int kk;
	for (ii=0; ii < tab->n_r; ii++){
		// TODO: SHOULD WE ONLY INTERPOLATE ONLY THE VALUES WE NEED??? //
		// only interpolate values where we need them (radius is defined invers!)
		if (ii<=ind_rmin || ii>=ind_rmax+1){
			interpol_a_mu0(ii, ifac_a, ifac_mu0, ind_a, ind_mu0, cached_tab_sysPar,tab);
		} else {  // set everything we won't need to 0 (just to be sure)
			cached_tab_sysPar->gmin[ii]=0.0;
			cached_tab_sysPar->gmax[ii]=0.0;
		    for (jj=0; jj<tab->n_g;jj++){
			    for (kk=0; kk<2;kk++){
			    	cached_tab_sysPar->trff[ii][jj][kk]=0.0;
			    	cached_tab_sysPar->cosne[ii][jj][kk]=0.0;
			    }
			 }
		}
	}

	/****************************/
	/** 2 **  Bin to Fine Grid **/
	/****************************/

	relSysPar* sysPar = (*sysPar_inp);
	// only need to initialize and allocat memory if not already loaded
	if (sysPar == NULL){
		sysPar = new_relSysPar(N_FRAD,tab->n_g,status);
		CHECK_STATUS_VOID(*status);
	}
	get_fine_radial_grid(rin,rout,sysPar);

	// let's try to be as efficient as possible here (note that "r" DEcreases)
	assert(ind_rmin>0); // as defined inverse, re[ind_rmin+1] is the lowest value
	assert((cached_tab_sysPar->re[ind_rmin+1]<=rin));
	assert((cached_tab_sysPar->re[ind_rmin]>=rin));
	assert((cached_tab_sysPar->re[ind_rmax+1]<=rout));
	assert((cached_tab_sysPar->re[ind_rmax]>=rout));
	assert(ind_rmax <= ind_rmin);
	assert(rout<=1000.0);

	double ifac_r;
	int ind_tabr=ind_rmin;

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
		// assert(ifac_r>=0.0);

		// we only allow extrapolation (i.e. ifac_r < 0) for the last bin
		if (ifac_r >1.0 && ind_tabr>0){
			RELXILL_ERROR("interpolation of rel_table on fine radial grid failed due to corrupted grid",status);
			printf("   --> radius %.4e not found in [%.4e,%.4e]  \n",
					sysPar->re[ii],cached_tab_sysPar->re[ind_tabr+1],cached_tab_sysPar->re[ind_tabr]);
			CHECK_STATUS_VOID(*status);
		}

		int jj; int kk;
		for (jj=0; jj<sysPar->ng; jj++){
			for (kk=0; kk<2; kk++){
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
	(*sysPar_inp) = sysPar;

	return;
}



/* function to get the system parameters; decides if values need to be re-computed
 * interpolate values loaded from the table if necessary */
static relSysPar* get_system_parameters(relParam* param, int* status){

	// only re-do the interpolation if rmin,rmax,a,mu0 changed
	// or if the cached parameters are NULL
	if (redo_get_system_parameters(param,cached_rel_param)){
		double mu0 = cos(param->incl);
		interpol_relTable(&cached_relSysPar,param->a,mu0,param->rin,param->rout,status);
		CHECK_STATUS_RET(*status,NULL);

		// get emissivity profile
		calc_emis_profile(param, cached_rel_param,cached_relSysPar, status);
		CHECK_STATUS_RET(*status,NULL);
	}

	return cached_relSysPar;
}

/** get new structure to store the relline spectrum (possibly for several zones)
    important note: ener has n_ener+1 number of bins **/
rel_spec* new_rel_spec(int nzones, const int n_ener, int*status){

	rel_spec* spec = (rel_spec*) malloc(sizeof(rel_spec));
	CHECK_MALLOC_RET_STATUS(spec,status,NULL);


	spec->n_zones = nzones;
	spec->n_ener  = n_ener;

	spec->flux = (double**) malloc (spec->n_zones * sizeof(double*) );
	CHECK_MALLOC_RET_STATUS(spec->flux,status,spec);

	int ii;


	for (ii=0; ii<spec->n_zones; ii++){
		spec->flux[ii] = (double*) malloc ( n_ener * sizeof(double) );
		CHECK_MALLOC_RET_STATUS(spec->flux[ii],status,spec);
	}

	spec->rgrid = (double*) malloc ( (spec->n_zones+1) * sizeof(double) );
	CHECK_MALLOC_RET_STATUS(spec->rgrid,status,spec);

	spec->ener = (double*) malloc ( (spec->n_ener+1) * sizeof(double) );
	CHECK_MALLOC_RET_STATUS(spec->ener,status,spec);


	spec->rel_cosne = NULL;


	return spec;
}

/** get new structure to store the cosne distribution spectrum (possibly for several zones) **/
rel_cosne* new_rel_cosne(int nzones, int n_incl, int*status){

	rel_cosne* spec = (rel_cosne*) malloc(sizeof(rel_cosne));
	CHECK_MALLOC_RET_STATUS(spec,status,NULL);


	spec->n_zones = nzones;
	spec->n_cosne = n_incl;

	spec->cosne = (double*) malloc (spec->n_cosne * sizeof(double) );
	CHECK_MALLOC_RET_STATUS(spec->cosne,status,spec);


	spec->dist = (double**) malloc (spec->n_zones * sizeof(double*) );
	CHECK_MALLOC_RET_STATUS(spec->dist,status,spec);

	int ii;

	for (ii=0; ii<spec->n_zones; ii++){
		spec->dist[ii] = (double*) malloc ( n_incl * sizeof(double) );
		CHECK_MALLOC_RET_STATUS(spec->dist[ii],status,spec);
	}

	return spec;
}



/** initialize the rel_spec structure **/
static void init_rel_spec(rel_spec** spec, relParam* param, xillTable* xill_tab,
		double** pt_ener, const int n_ener, int* status ){

	/** in case of the relxill-LP model multiple zones are used **/
	int nzones = get_num_zones(param->model_type, param->emis_type);

	if ((*spec)==NULL){
		(*spec) = new_rel_spec(nzones,n_ener,status);
	} else {
		// check if the number of zones changed or number of energy bins
		if ( (nzones != (*spec)->n_zones ) || (n_ener != (*spec)->n_ener )){
		// -> if yes, we need to free memory and re-allocate the appropriate space
		free_rel_spec( (*spec) );
		(*spec) = new_rel_spec(nzones,n_ener,status);
		}
	}

	int ii;
	for (ii=0; ii<=(*spec)->n_ener; ii++ ){
		(*spec)->ener[ii] = (*pt_ener)[ii];
	}


	if (xill_tab!=NULL){
		if ((*spec)->rel_cosne == NULL) {
			(*spec)->rel_cosne = new_rel_cosne(nzones,xill_tab->n_incl,status);
		}
		int ii;
		for (ii=0; ii<(*spec)->rel_cosne->n_cosne; ii++ ){
			(*spec)->rel_cosne->cosne[ii] = cos(xill_tab->incl[ii]*M_PI/180);
		}
	}

	get_rzone_grid(param->rin, param->rout, (*spec)->rgrid, nzones, status);
	CHECK_STATUS_VOID(*status);

	return;
}

static void zero_rel_spec_flux(rel_spec* spec){
	int ii; int jj;
	for (ii=0; ii<spec->n_zones;ii++){
		for (jj=0; jj<spec->n_ener;jj++){
			spec->flux[ii][jj] = 0.0;
		}
		if ( spec->rel_cosne!=NULL ){
			for (jj=0; jj<spec->rel_cosne->n_cosne;jj++){
				spec->rel_cosne->dist[ii][jj]=0.0;
			}
		}
	}
}


/** relat. transfer function, which we will need to integrate over the energy bin then **/
static str_relb_func* new_str_relb_func(relSysPar* sysPar, int* status){
	str_relb_func* str = (str_relb_func*) malloc(sizeof(str_relb_func));
	CHECK_MALLOC_RET_STATUS(str,status,NULL);

	str->gstar = sysPar->gstar;
	str->ng = sysPar->ng;

	//str->cache_val_relb_func = (double*) malloc(2*sizeof(double));
	//CHECK_MALLOC_RET_STATUS(str->cache_val_relb_func,status,str);

	return str;
}

/** relat. function which we want to integrate **/
static double relb_func(double eg, int k, str_relb_func* str){

  // get the redshift from the energy
  // double eg = e/line_energy;
  double egstar = (eg-str->gmin)*str->del_g;

  // find the indices in the original g-grid, but check first if they have already been calculated
  int ind;
  if (!((egstar>=str->gstar[str->save_g_ind])&&(egstar<str->gstar[str->save_g_ind+1]))){
	  str->save_g_ind  = binary_search(str->gstar,str->ng,egstar);
  }
  ind = str->save_g_ind;

  double inte = (egstar - str->gstar[ind])/(str->gstar[ind+1] - str->gstar[ind]);
  double inte1=1.0-inte;
  double ftrf = inte*str->trff[ind][k] + inte1*str->trff[ind+1][k];

  // double fmu0 = inte*str->cosne[ind,k] + inte1*str->cosne[ind+1,k];

  /** isotropic limb law by default (see Svoboda (2009))
  limb = 1.D0
  if (limb_law.eq.1) then   !Laor(1991)
       limb = (1.D0 + 2.06D0*fmu0)
  else if (limb_law.eq.2) then  !Haardt (1993)
     limb = log(1.D0+1.D0/fmu0)
  endif **/

 // ! INTEGRAL !
 return pow((eg*str->re),3)/((str->gmax-str->gmin)*sqrt(egstar - egstar*egstar))*ftrf*str->emis;
}

/** Romberg Integration Routine **/
static double RombergIntegral(double a,double b,int k, str_relb_func* str, double line_ener){
  const double prec = 0.05;
  double obtprec = 1.0;
  const int itermin = 0;
  int itermax = 5;
  const int maxiter = 6;
  double t[maxiter+1][maxiter+1];

  int niter = 0;

  if (itermax>maxiter) {
	  itermax=maxiter;
  }

  // check if this value has already been calculated
  double r;
  if (str->cached_relbf) {
     r = str->cache_val_relb_func[k];
  } else {
     r = relb_func(a,k,str);
  }

  str->cache_val_relb_func[k] = relb_func(b,k,str);
  str->cache_rad_relb_fun = str->re;
  // rb(k) = RELB_FUNC(b,k);


  int ii;

  double ta = (r + str->cache_val_relb_func[k]) / 2.0;
  double pas=b-a;
  double pasm=1.0;
  double s=ta;
  t[0][0]=ta*pas;
  while ( (niter<itermin) || ((obtprec > prec) && (niter <= itermax) )) {
	  niter++;
	  pas=pas/2.0;
	  pasm=pasm/2.0;
	  s=ta;
	  for (ii=1; ii<=pow(2,niter)-1; ii++) {
		  s += relb_func(a+pas*ii,k,str);
	  }
	  t[0][niter]=s*pas;
	  r=1.0;
	  for (ii=1; ii<=niter; ii++){
		  r *= 4.0;
		  int jj=niter-ii;
		  t[ii][jj] = (r*t[ii-1][jj+1] - t[ii-1][jj])/(r-1.0);
	  }
	  obtprec = fabs(t[niter][0] - t[niter-1][0])/t[niter][0];
 }

  return t[niter][0];
}


static void free_str_relb_func(str_relb_func* str){
	if (str!=NULL){
//		free(str->cache_val_relb_func);
		free(str);
	}
}


/** function which makes an approximated integration for gstar->0/1
    this is only done within gstar=[0,H] and gstar[H,1-H]
    input:   bin_lo and bin_hi
    output:  area of the bin (= luminosity = E/dt/bin) )  **/
static double int_edge(double blo, double bhi, double h, str_relb_func* str, double line_energy){


  // get the value of the Luminosity on the point closest to the ones to be approximated ("H")

	double hex;
	double lo; double hi;

 // #1: upper or lower limit -> write the corresponding gstar-value to hex
  if (blo <= 0.5) {
    hex = h;
    lo = blo;
    hi = bhi;
  } else {
    hex = 1.0-h;
    /**  variable transformation for upper limit x -> 1 - x
         leads to the same integral
         (if the correct normalization is chosen; hex keeps track of that) **/
    lo = 1.0-bhi;
    hi = 1.0-blo;
  }

  // #2: get the normalization value
  int k=0;
  double norm = 0.0;
  for (k=0; k<2; k++) {
	  norm = norm + relb_func(gstar2ener(hex,str->gmin,str->gmax,line_energy),k,str);
  }


 // #3: thanks to variable transformation:
 //   both cases are described by the one for [0,H]
  norm = norm*sqrt(h);

  return 2*norm*(sqrt(hi)-sqrt(lo))*line_energy*(str->gmax-str->gmin);
}


/** function which calculates the normal integral by romberg's method
(it is acutally jsut a wrapper which sets the correct parameters
and coordinates the integration over k=1,2)
input:   bin_lo and bin_hi
output:  area of the bin (= luminosity = E/dt/bin) ) **/
static double int_romb(double lo, double hi, str_relb_func* str, double line_energy){

	double flu = 0.0;

	/** Problems with the blue peak WITHOUT Romberg --> Romberg for blue peak only **/
	int k;
	if (lo>=line_energy){
		for (k=0;k<2;k++){
			flu += RombergIntegral(lo,hi,k,str,line_energy);
		}
	} else {
		for (k=0;k<2;k++){
			flu += relb_func((hi+lo)/2.0,k,str)*(hi-lo);
		}
	}

	return flu;
}

/** integrate the flux bin (see Dauser+2010, MNRAS for details) **/
static double integ_relline_bin(str_relb_func* str, double rlo0, double rhi0){

	double line_ener = 1.0;
	double flu=0.0;

	double gblo = (rlo0/line_ener - str->gmin)*str->del_g;
	if (gblo<0.0) {
		gblo=0.0;
	} else if (gblo>1.0){
		gblo = 1.0;
	}

	double gbhi = (rhi0/line_ener - str->gmin)*str->del_g;
	if (gbhi<0.0) {
		gbhi=0.0;
	} else if (gbhi>1.0){
		gbhi = 1.0;
	}
	if (gbhi==0){
		return 0.0;
	}

	double rlo = rlo0;
	double rhi = rhi0;

	double hlo; double hhi;
	// #1: low approx. integration
	if (gblo<=GFAC_H) {
		// range of the low-app.-integration
		hlo = gblo;
		hhi = GFAC_H;
		// begin of the 'real' integration
		rlo = gstar2ener(GFAC_H,str->gmin,str->gmax,line_ener);
		// .. but also check if this integration is necessary
		if (gbhi<=GFAC_H) {
			hhi = gbhi;
			rlo = -1.0;
		}
		// approximated integration
		flu = flu + int_edge(hlo,hhi,GFAC_H,str,line_ener);
	}

	// #2: upper limit approximation to be taken into account?
	if (gbhi >= (1.0-GFAC_H)) {
		// range of the upper-app.-integration
		hhi = gbhi;
		hlo = 1.0-GFAC_H;

		// begin of the 'real' integration
		rhi = gstar2ener(1-GFAC_H,str->gmin,str->gmax,line_ener);
		// .. but also check if this integration is necessary
		if (gblo>=(1.0-GFAC_H)) {
			hlo = gblo;
			rhi = -1.0;
		}

		/** the same approximated integration as in case #1 can be
	           applied, if one makes a variable transformation **/
		flu = flu + int_edge(hlo,hhi,GFAC_H,str,line_ener);
	}

	// #3: real integration (only if necessary)
	if ((rhi>=0) && (rlo>=0)) {

		// has the function relb_func been calculated at the lower bin boundary before?
		// (should've been upper bound before; also make sure we haven't changed the radial bin!)
		if ( (fabs(rlo - str->cache_bin_ener) < cache_limit) && (fabs(str->re - str->cache_rad_relb_fun) < cache_limit)) {
			str->cached_relbf = 1;
		} else {
			str->cached_relbf = 0;
		}
		flu = flu  + int_romb(rlo,rhi,str,line_ener);
	}


	return flu;
}

static void	set_str_relbf(str_relb_func* str, double re, double gmin, double gmax, double** trff, double** cosne, double emis){
	str->re = re;
	str->gmin = gmin;
	str->gmax = gmax;
	str->del_g = 1./(gmax-gmin);
	str->emis = emis;

	str->trff = trff;
	str->cosne = cosne;

	str->cache_bin_ener = -1.0;
	str->cached_relbf = 0;

	str->save_g_ind = 0;
}


/** function to properly re-normalize the relline_profile **/
static void renorm_relline_profile(rel_spec* spec, relParam* rel_param){
	// normalize to 'cts/bin'
	int ii; int jj;
	double sum = 0.0;
	for (ii=0; ii<spec->n_zones; ii++){
		for (jj=0; jj<spec->n_ener; jj++){
			spec->flux[ii][jj] /= 0.5*( spec->ener[jj] + spec->ener[jj+1] );
			sum += spec->flux[ii][jj];
		}
	}

	/** only renormalize if not the relxill model **/
	if ((! is_relxill_model(rel_param->model_type) ) || (rel_param->emis_type != EMIS_TYPE_LP)) {
		for (ii=0; ii<spec->n_zones; ii++){
			for (jj=0; jj<spec->n_ener; jj++){
				spec->flux[ii][jj] /= sum;
			}
		}
	}

	if (spec->rel_cosne!=NULL){
		for (ii=0; ii<spec->n_zones; ii++){
			// normalize it for each zone, the overall flux will be taken care of by the normal structure
			sum = 0.0;
			for (jj=0; jj<spec->rel_cosne->n_cosne; jj++){
				sum += spec->rel_cosne->dist[ii][jj];
			}
			for (jj=0; jj<spec->rel_cosne->n_cosne; jj++){
				spec->rel_cosne->dist[ii][jj] /= sum;
			}
		}
	}
}

/** get the bin for a certain emission angle (between [0,n_incl-1] **/
int static get_cosne_bin(double mu, rel_cosne* dat){
	return (( int ) (dat->n_cosne*(1 - mu) + 1)) - 1 ;
}

/** calculate the relline profile(s) for all given zones **/
str_relb_func* cached_str_relb_func = NULL;
void relline_profile(rel_spec* spec, relSysPar* sysPar, int* status){

	double line_ener=1.0;

	// very important: set all fluxes to zero
	zero_rel_spec_flux(spec);

	if (cached_str_relb_func==NULL){
		cached_str_relb_func = new_str_relb_func(sysPar, status);
	}


	int ii; int jj;
	for (ii=0; ii<sysPar->nr; ii++){
	     // gstar in [0,1] + corresponding energies (see gstar2ener for full formula)
	     double egmin = sysPar->gmin[ii]*line_ener;
	     double egmax = sysPar->gmax[ii]*line_ener;

	     // check if the expected energy-bins are needed
	     if ((egmax>spec->ener[0]) && (egmin<spec->ener[spec->n_ener])){


	        /**  make sure that integration is only done inside the
	             given energy range **/
	        if (egmin<spec->ener[0]){
	           egmin = spec->ener[0];
	        }
	        if (egmax>spec->ener[spec->n_ener]){
	           egmax=spec->ener[spec->n_ener];
	        }

	        /** search for the indices in the ener-array
	            index is such that: ener[k]<=e<ener[k+1] **/
	        int ielo = binary_search(spec->ener,spec->n_ener+1,egmin);
	        int iehi = binary_search(spec->ener,spec->n_ener+1,egmax);

	        // in which ionization bin are we?
	        int izone = binary_search(spec->rgrid,spec->n_zones+1,sysPar->re[ii]);

	        // trapez integration
	       	double weight = trapez_integ_single(sysPar->re,ii,sysPar->nr);

	       	// set the current parameters in a cached structure (and reset some values) [optimizes speed]
	       	set_str_relbf(cached_str_relb_func,
	       			sysPar->re[ii], sysPar->gmin[ii],sysPar->gmax[ii],
	       			sysPar->trff[ii],sysPar->cosne[ii],
					sysPar->emis[ii]);

	       	// lastly, loop over the energies
	       	for (jj=ielo; jj<=iehi; jj++){
	       		double tmp_var = integ_relline_bin(cached_str_relb_func,spec->ener[jj],spec->ener[jj+1]);
	       		spec->flux[izone][jj] += tmp_var*weight;
	       	}

	       	/** only calculate the distribution if we need it here  **/
	       	if (spec->rel_cosne != NULL){
	       		int kk; int imu;
	       		str_relb_func* da = cached_str_relb_func; // define a shortcut
	       		for (jj=0; jj<sysPar->ng; jj++){
	       			double g = da->gstar[jj]*(da->gmax-da->gmin) + da->gmin;
	       			for (kk=0; kk<2; kk++){
	       				imu = get_cosne_bin(da->cosne[jj][kk],spec->rel_cosne);

	       				spec->rel_cosne->dist[izone][imu] +=
	       						da->re*pow(2*M_PI*g*da->re,2)/
	       		                sqrt(da->gstar[jj] - da->gstar[jj]*da->gstar[jj])*
								da->trff[jj][kk]* da->emis
								* weight * sysPar->d_gstar[jj] ;
	       			}
	       		}
	       	} /** end calculating angular distribution **/


		}
	}

}

/** convolve the (bin-integrated) spectra f1 and f2 (which need to have a certain binning)
 *  fout: gives the output
 *  f1 input (reflection) specrum
 *  f2 filter
 *  ener has length n+1 and is the energy array
 * **/
void fft_conv_spectrum(double* ener, double* f1, double* f2, double* fout, int n, int* status){

	long m=0;
	switch(n)	{
		case 512:  m = 9; break;
		case 1024: m = 10; break;
		case 2048: m = 11; break;
		case 4096: m = 12; break;
		case 8192: m = 13; break;
		default: *status=EXIT_FAILURE; printf(" *** error: Number of Bins %i not allowed in Convolution!! \n",n);break;
	}
	CHECK_STATUS_VOID(*status);

	/* need to find out where the 1keV for the filter is, which defines if energies are blue or redshifted*/
	if (save_1eV_pos==0 ||
		(! ( (ener[save_1eV_pos]<=1.0) &&
		(ener[save_1eV_pos+1]>1.0)  ))){
		save_1eV_pos = binary_search(ener,n+1,1.0);
	}


	int ii;
	int irot;
	double x1[n]; double y1[n];
	double x2[n]; double y2[n];
	double xcomb[n]; double ycomb[n];
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum_fft = 0.0;
	for (ii=0; ii<n; ii++){
		x1[ii] = f1[ii]/(ener[ii+1]-ener[ii]);
		irot = (ii-save_1eV_pos+n) % n ;
		x2[irot ] = f2[ii]/(ener[ii+1]-ener[ii]);
		y1[ii] = 0.0;
		y2[ii] = 0.0;
		sum2 += f2[ii];
		sum1 += f1[ii];
	}

	/** TODO: we could cache either the relat. or the xillver part, as only one of the two changes (reduce time by 1/3 for convolution) **/
	FFT_R2CT(1, m, x1, y1);
	FFT_R2CT(1, m, x2, y2);

	/* complex multiplication
	 * (we need the real part, so we already use the output variable here
	 *  to save computing time */
	for (ii=0; ii<n; ii++){
		xcomb[ii] = x1[ii]*x2[ii] - y1[ii]*y2[ii];
		ycomb[ii] = y1[ii]*x2[ii] + x1[ii]*y2[ii];
	}

	FFT_R2CT(-1, m, xcomb, ycomb);

	for (ii=0; ii<n; ii++){
		fout[ii] = xcomb[ii] * (ener[ii+1]-ener[ii]);
		sum_fft += fout[ii];
	}

	for (ii=0; ii<n; ii++){
//		fout[ii] *= sum1*sum2 / sum_fft;
	}
}


static void print_angle_dist(rel_cosne* spec, int izone){


	FILE* fp =  fopen ( "testfiles/test_angle_dist.dat","w+" );
	int ii;
	for(ii=0; ii<spec->n_cosne; ii++ ){
		fprintf(fp," mu=%e %e \n",spec->cosne[ii],
				// acos(spec->cosne[ii])*180/M_PI,
				spec->dist[izone][ii]);
	}
	if (fclose(fp)) exit(1);

}


static void calc_xillver_angdep(double* xill_flux, xill_spec* xill_spec,
		double* cosne, double* dist, int n_incl, int* status){

	int ii; int jj;
	for (ii=0; ii<xill_spec->n_ener;ii++){
		xill_flux[ii] = 0.0;
	}

	double mufac;
	for (ii=0; ii<xill_spec->n_incl;ii++){
		/**  0.5*F*mu*dmu (Javier Note, Eq. 25) **/
		mufac = dist[ii];   // actually it is also mutliplied by  dmu*nincl = 1/nincl * nincl = 1
		for (jj=0; jj<xill_spec->n_ener;jj++){
			xill_flux[jj] += mufac*xill_spec->flu[ii][jj];
		}
	}

}

static double get_rzone_energ_shift(relParam* param, double rad){

	// can be implemented by using the del_emit from relSysParam in principle
	if (param->beta){
		printf(" *** error :  energy shift of the incident spectrum for beta>0 not implemented yet ; assuming beta=0 here ");
	}

	// as beta=0, we can choose any value for del_emit

	return gi_potential_lp(rad,param->a,param->height,0.0,M_PI/2);
}

/** BASIC RELCONV FUNCTION : convole any input spectrum with the relbase kernel
 *  (ener has the length n_ener+1)
 *  **/
void relconv_kernel(double* ener_inp, double* spec_inp, int n_ener_inp, relParam* rel_param, int* status ){

	/* get the (fixed!) energy grid for a RELLINE for a convolution
	 * -> as we do a simple FFT, we can now take into account that we
	 *    need it to be number = 2^N */

	// always do the convolution on this grid
	/** only do the calculation once **/
	if (ener_std == NULL){
		ener_std = (double*) malloc( (n_ener_std+1) * sizeof(double));
		CHECK_MALLOC_VOID_STATUS(ener_std,status);
		get_log_grid(ener_std, (n_ener_std+1), EMIN_RELXILL, EMAX_RELXILL);
	}
	int n_ener = n_ener_std;
	double* ener = ener_std;

	rel_spec* rel_profile = relbase(ener, n_ener, rel_param, NULL, status);

	// simple convolution only makes sense for 1 zone !
	assert(rel_profile->n_zones==1);

	double rebin_flux[n_ener];
	double conv_out[n_ener];
	rebin_spectrum(ener,rebin_flux,n_ener,
			ener_inp, spec_inp, n_ener_inp );

	// convolve the spectrum
	fft_conv_spectrum(ener, rebin_flux, rel_profile->flux[0], conv_out,  n_ener, status);
	CHECK_STATUS_VOID(*status);

	// rebin to the output grid
	rebin_spectrum(ener_inp,spec_inp,n_ener_inp, ener, conv_out, n_ener);

}


/** BASIC RELXILL KERNEL FUNCTION : convole a xillver spectrum with the relbase kernel
 * (ener has the length n_ener+1)
 *  **/

void relxill_kernel(double* ener_inp, double* spec_inp, int n_ener_inp, xillParam* xill_param, relParam* rel_param, int* status ){

	/* get the (fixed!) energy grid for a RELLINE for a convolution
	 * -> as we do a simple FFT, we can now take into account that we
	 *    need it to be number = 2^N */


	/** only do the calculation once **/
	if (ener_std == NULL){
		ener_std = (double*) malloc( (n_ener_std+1) * sizeof(double));
		CHECK_MALLOC_VOID_STATUS(ener_std,status);
		get_log_grid(ener_std, (n_ener_std+1), EMIN_RELXILL, EMAX_RELXILL);
	}
	int n_ener = n_ener_std;
	double* ener = ener_std;


	xillTable* xill_tab = NULL;
	get_init_xillver_table(&xill_tab,xill_param,status);

	xill_spec* xill_spec = NULL;

	// todo: can/should we add caching here directly??
	rel_spec* rel_profile = relbase(ener, n_ener, rel_param, xill_tab, status);

	if (is_debug_run()){
		save_xillver_spectrum(ener,rel_profile->flux[0],n_ener,"testfiles/test_relxill_conv_spectrum.dat");
	}


	// make sure the output array is set to 0
	int ii;
	for (ii=0; ii<n_ener_inp;ii++){
		spec_inp[ii]=0.0;
	}

	/*** LOOP OVER THE RADIAL ZONES ***/
	double xill_flux[n_ener];
	double xill_angdist_inp[xill_tab->n_ener];
	double conv_out[n_ener];
	double single_spec_inp[n_ener_inp];
	for (ii=0; ii<n_ener;ii++){
		conv_out[ii] = 0.0;
	}

	double ecut0 = xill_param->ect;

	int jj;
	for (ii=0; ii<rel_profile->n_zones;ii++){

		double test_sum_relline = 0.0;
		double test_sum_relxill = 0.0;
		double test_sum_xillver = 0.0;

		// now calculate the reflection spectra for each zone (using the angular distribution)
		assert(rel_profile->rel_cosne != NULL);
		if (is_debug_run()==1){
			print_angle_dist(rel_profile->rel_cosne,ii);
		}

		// get the energy shift in order to calculate the proper Ecut value (if nzones>1)
		// (the latter part of the IF is a trick to get the same effect as NZONES=1 if during a running
		//  session the number ofd zones is reset)
		if (rel_profile->n_zones==1 || get_num_zones(rel_param->model_type,rel_param->emis_type)==1){
			xill_param->ect = ecut0;
		} else {
			// choose the (linear) middle of the radial zone
			double rzone = 0.5*(rel_profile->rgrid[ii]+rel_profile->rgrid[ii+1]);
			xill_param->ect = ecut0 * ( 1 + grav_redshift(rel_param) ) * get_rzone_energ_shift(rel_param,rzone ) ;
			// printf(" rzone=%.2e-%.2e -> Ecut = %.3f (ecut0=%.2e)\n",rel_profile->rgrid[ii],rel_profile->rgrid[ii+1],xill_param->ect,ecut0);
			// printf("  --> energ_shift %.5e  --  grav redshift %.5e \n",get_rzone_energ_shift(rel_param,rzone ),  ( 1 + grav_redshift(rel_param) ));
		}

		// call the function which calculates the xillver spectrum
		xill_spec = get_xillver_spectra(xill_param,status);
		CHECK_STATUS_VOID(*status);

		// get angular distribution
		calc_xillver_angdep(xill_angdist_inp,xill_spec,rel_profile->rel_cosne->cosne,
				rel_profile->rel_cosne->dist[ii],rel_profile->rel_cosne->n_cosne,status);

		rebin_spectrum(ener,xill_flux,n_ener,
				xill_spec->ener, xill_angdist_inp, xill_spec->n_ener );

		// convolve the spectrum
		fft_conv_spectrum(ener, xill_flux, rel_profile->flux[ii], conv_out,  n_ener, status);
		CHECK_STATUS_VOID(*status);

		// rebin to the output grid
		rebin_spectrum(ener_inp,single_spec_inp,n_ener_inp, ener, conv_out, n_ener);


		for (jj=0; jj<n_ener;jj++){
			if (ener[jj] > EMIN_XILLVER && ener[jj+1] < EMAX_XILLVER  ){
				test_sum_relline += rel_profile->flux[ii][jj];
				test_sum_relxill += conv_out[jj];
				test_sum_xillver += xill_flux[jj];
			}
		}

		/** avoid problems where no relxill bin falls into an ionization bin **/
		if (test_sum_relline < 1e-12){
			continue;
		}

		// add it to the final output spectrum
		for (jj=0; jj<n_ener_inp;jj++){
			spec_inp[jj] += single_spec_inp[jj]*test_sum_relline*test_sum_xillver/test_sum_relxill;
		}

		if (is_debug_run()){
			char* vstr;
			double test_flu[n_ener_inp];
			for (jj=0; jj<n_ener_inp;jj++){
				test_flu[jj] = single_spec_inp[jj]*test_sum_relline*test_sum_xillver/test_sum_relxill;
			}
			if (asprintf(&vstr, "testfiles/test_relxill_spec_zones_%03i.dat", ii+1) == -1){
				RELXILL_ERROR("failed to get filename",status);
			}
			save_xillver_spectrum(ener_inp,test_flu,n_ener_inp,vstr);
		}

	}


	// lastely, we make the spectrum normalization independent of the ionization parameter
	// todo: put this in a ionization gradient routine ??
	for (ii=0; ii<n_ener_inp; ii++){
		spec_inp[ii] /= pow(10,xill_param->lxi);
		if (fabs(xill_param->dens - 15) > 1e-6 ){
			spec_inp[ii] /= pow(10,xill_param->dens - 15);
		}
	}


	/** add a primary spectral component and normalize according to the given refl_frac parameter**/
	add_primary_component(ener_inp,n_ener_inp,spec_inp,rel_param,xill_param, status);

	free_xill_spec(xill_spec);

	return;
}

/** function adding a primary component with the proper norm to the flux **/
void add_primary_component(double* ener, int n_ener, double* flu, relParam* rel_param,
		xillParam* xill_param, int* status){

	double pl_flux[n_ener];

	int ii;
	double pl_flux_xill[n_ener_xill];

	/** need to create an energy grid for the primary component to fulfill the XILLVER NORM condition (Dauser+2016) **/
	if (ener_xill == NULL){
		ener_xill = (double*) malloc( (n_ener_xill+1) * sizeof(double));
		CHECK_MALLOC_VOID_STATUS(ener_xill,status);
		get_log_grid(ener_xill,n_ener_xill+1,ener_xill_norm_lo,ener_xill_norm_hi);
	}

	/** 1 **  calculate the primary spectrum  **/
	if (xill_param->prim_type == PRIM_SPEC_ECUT){

		// TODO: change to new defintion ???
		double ecut_rest = xill_param->ect;

		for (ii=0; ii<n_ener_xill; ii++){
			pl_flux_xill[ii] = exp(1.0/ecut_rest) *
		             pow(0.5*(ener_xill[ii]+ener_xill[ii+1]),-xill_param->gam) *
		             exp( -0.5*(ener_xill[ii]+ener_xill[ii+1])/ecut_rest) *
		             (ener_xill[ii+1] - ener_xill[ii]);
		}

	} else if ( xill_param->prim_type == PRIM_SPEC_NTHCOMP){
		printf(" **** warning: NTHCOMP primary continuum not yet implemented \n");
	} else {
		RELXILL_ERROR("trying to add a primary continuum to a model where this does not make sense (should not happen!)",status);
		return;
	}



	/** 2 **  get the normalization of the spectrum with respect to xillver **/
	double norm_xill = pow(10,xill_param->dens) / 4.0 / M_PI;
	double keV2erg = 1.602177e-09;

	double sum_pl = 0.0;
	for (ii=0; ii<n_ener_xill; ii++){
	     sum_pl += pl_flux_xill[ii] * 0.5*(ener_xill[ii] + ener_xill[ii+1]) * 1e20 * keV2erg;
	}
	double norm_pl = norm_xill / sum_pl;   // normalization defined, e.g., in Appendix of Dauser+2016

	/** bin the primary continuum onto the given grid **/
	rebin_spectrum( ener, pl_flux,n_ener, ener_xill, pl_flux_xill, n_ener_xill);

	/** 2 **  decide if we need to do relat. calculations **/
	if (is_xill_model(xill_param->model_type) ){

		for (ii=0; ii<n_ener; ii++){
			pl_flux[ii] *= norm_pl;
			flu[ii] *= fabs(xill_param->refl_frac);
		}
	} else {

		assert(rel_param!=NULL);

		// should be cached, as it has been calculated before (todo: check this!!)
		relSysPar* sysPar = get_system_parameters(rel_param, status);

		/** 3 **  calculate predicted reflection fraction and check if we want to use this value **/
		lpReflFrac* struct_refl_frac = calc_refl_frac(sysPar, rel_param,status);
		CHECK_STATUS_VOID(*status);

		if ((xill_param->fixReflFrac==1)||(xill_param->fixReflFrac==2)) {
			xill_param->refl_frac = struct_refl_frac->refl_frac;
		}

		/** 4 ** and apply it to primary and reflected spectra **/
		if (rel_param->emis_type == EMIS_TYPE_LP) {
			double g_inf = sqrt( 1.0 - ( 2*rel_param->height /
					(rel_param->height*rel_param->height + rel_param->a*rel_param->a)) );
			for (ii=0; ii<n_ener; ii++) {
				pl_flux[ii] *= norm_pl * pow(g_inf,xill_param->gam) * (struct_refl_frac->f_inf / 0.5);
				flu[ii] *= fabs(xill_param->refl_frac) / struct_refl_frac->refl_frac;
			}

		} else {
			for (ii=0; ii<n_ener; ii++){
				pl_flux[ii] *= norm_pl;
				flu[ii] *= fabs(xill_param->refl_frac);
			}
		}


		/** 5 ** if desired, we ouput the reflection fraction and strength (as defined in Dauser+2016) **/
		if ((xill_param->fixReflFrac == 2) && (rel_param->emis_type==EMIS_TYPE_LP)) {

			/** the reflection strength is calculated between RSTRENGTH_EMIN and RSTRENGTH_EMAX **/
			// todo: all this to be set by a qualifier

			int imin = binary_search(ener,n_ener+1,RSTRENGTH_EMIN);
			int imax = binary_search(ener,n_ener+1,RSTRENGTH_EMAX);

			sum_pl = 0.0;
			double sum = 0.0;
			for (ii=imin; ii<=imax; ii++){
				sum_pl += pl_flux[ii];
				sum += flu[ii];
			}

			printf(" the reflection fraction for a = %.2f and %.2f rg is: %.3f and the reflection strength is: %.3f \n",
					rel_param->a,rel_param->height,struct_refl_frac->refl_frac,sum/sum_pl);
		}

		/** free the reflection fraction structure **/
		free(struct_refl_frac);

	}




	/** 6 ** add power law component only if desired (i.e., refl_frac > 0)**/
	  if (xill_param->refl_frac >= 0) {
	     for (ii=0; ii<n_ener; ii++) {
	        flu[ii] += pl_flux[ii];
	     }
	  }

}

/** print the relline profile   **/
void save_relline_profile(rel_spec* spec){

	FILE* fp =  fopen ( "test_relline_profile.dat","w+" );
	int ii;
	for (ii=0; ii<spec->n_ener; ii++){
		fprintf(fp, " %e \t %e \t %e \n",spec->ener[ii],spec->ener[ii+1],spec->flux[0][ii]);
	}
	if (fclose(fp)) exit(1);
}

/** print the relline profile   **/
static void save_emis_profile(double* rad, double* intens, int n_rad){

	FILE* fp =  fopen ( "test_emis_profile.dat","w+" );
	int ii;
	for (ii=0; ii<n_rad; ii++){
		fprintf(fp, " %e \t %e \n",rad[ii],intens[ii]);
	}
	if (fclose(fp)) exit(1);
}

static int comp_single_param_val(double val1, double val2){
	if (fabs(val1-val2) <= cache_limit){
		return 0;
	} else {
		return 1;
	}
}

static void set_cached_rel_param(relParam* par,int* status){

	if (cached_rel_param==NULL){
		cached_rel_param = (relParam*) malloc ( sizeof(relParam));
		CHECK_MALLOC_VOID_STATUS(cached_rel_param,status);
	}

	cached_rel_param->a = par->a;
	cached_rel_param->emis1 = par->emis1;
	cached_rel_param->emis2 = par->emis2;
	cached_rel_param->gamma = par->gamma;
	cached_rel_param->height = par->height;
	cached_rel_param->incl = par->incl;
	cached_rel_param->beta = par->beta;

	cached_rel_param->z = par->z;
	cached_rel_param->lineE = par->lineE;

	cached_rel_param->emis_type = par->emis_type;
	cached_rel_param->model_type = par->model_type;

	cached_rel_param->rbr = par->rbr;
	cached_rel_param->rin = par->rin;
	cached_rel_param->rout = par->rout;

}

static int comp_rel_param(relParam* cpar, relParam* par){
	if (comp_single_param_val(par->a,cpar->a)) return 1;
	if (comp_single_param_val(par->emis1,cpar->emis1)) return 1;
	if (comp_single_param_val(par->emis2,cpar->emis2)) return 1;
	if (comp_single_param_val(par->gamma,cpar->gamma)) return 1;
	if (comp_single_param_val(par->height,cpar->height)) return 1;
	if (comp_single_param_val(par->incl,cpar->incl)) return 1;
	if (comp_single_param_val(par->beta,cpar->beta)) return 1;

	/** need this for the current code to work; could be optimized **/
	if (comp_single_param_val(par->z,cpar->z)) return 1;
	if (comp_single_param_val(par->lineE,cpar->lineE)) return 1;


	if (comp_single_param_val( (double) par->emis_type, (double) cpar->emis_type)) return 1;
	if (comp_single_param_val( (double) par->model_type, (double) cpar->model_type)) return 1;

	if (comp_single_param_val(par->rbr,cpar->rbr)) return 1;
	if (comp_single_param_val(par->rin,cpar->rin)) return 1;
	if (comp_single_param_val(par->rout,cpar->rout)) return 1;

	return 0;
}

/* check if values, which need a re-computation of the relline profile, have changed */
int redo_relbase_calc(relParam* param, int* status){

	int redo = 1;

	if (cached_rel_param==NULL){
	} else {
		redo = comp_rel_param(cached_rel_param,param);
	}

	return redo;
}

/* the relbase function calculating the basic relativistic line shape for a given parameter setup
 * (assuming a 1keV line, by a grid given in keV!)
 * input: ener(n_ener), param
 * optinal input: xillver grid
 * output: photar(n_ener)     */
rel_spec* relbase(double* ener, const int n_ener, relParam* param, xillTable* xill_tab,int* status){

	// check caching here and also re-set the cached parameter values
	// TODO: also check if the energy grid changed!
	int redo = redo_relbase_calc(param,status);
	CHECK_STATUS_RET(*status,NULL);

	if (redo){

		// initialize parameter values
		relSysPar* sysPar = get_system_parameters(param,status);
		CHECK_STATUS_RET(*status,NULL);

		save_emis_profile(sysPar->re, sysPar->emis, sysPar->nr);

		// init the spectra where we store the flux
		init_rel_spec(&cached_rel_spec, param, xill_tab, &ener, n_ener, status);
		CHECK_STATUS_RET(*status,NULL);

		// calculate line profile (returned units are 'cts/bin')
		relline_profile(cached_rel_spec, sysPar, status);
		CHECK_STATUS_RET(*status,cached_rel_spec);

		// normalize it and calculate the angular disttribution (if necessary)
		renorm_relline_profile(cached_rel_spec,param);

		// store parameters such that we know what we calculated
		set_cached_rel_param(param,status);
		CHECK_STATUS_RET(*status,NULL);
	}

	// TODO: set photar, such that we can directly use this spectrum???

	assert(cached_rel_spec!=NULL);
	return cached_rel_spec;
}



void free_rel_cosne(rel_cosne* spec){
	if (spec!=NULL){
	//	free(spec->ener);  we do not need this, as only a pointer for ener is assigned
		free(spec->cosne);
		if (spec->dist!=NULL){
			int ii;
			for (ii=0; ii<spec->n_zones; ii++){
				free(spec->dist[ii]);
			}
		}
		free(spec->dist);
		free(spec);
	}
}

void free_rel_spec(rel_spec* spec){
	if (spec!=NULL){
		free(spec->ener);  // we do not need this, as only a pointer for ener is assigned
		free(spec->rgrid);
		if (spec->flux!=NULL){
			int ii;
			for (ii=0; ii<spec->n_zones; ii++){
				free(spec->flux[ii]);
			}
		}
		free(spec->flux);
		if (spec->rel_cosne != NULL){
			free_rel_cosne(spec->rel_cosne);
		}
		free(spec);
	}
}

void free_cached_tables(void){
	free_relTable(relline_table);
	free_relSysPar(cached_relSysPar);
	free_relSysPar(cached_tab_sysPar);
	free_cached_lpTable();

	free_cached_xillTable();

	free(cached_rel_param);

	free_rel_spec(cached_rel_spec);

	free_str_relb_func(cached_str_relb_func);

	free(ener_std);
	free(ener_xill);
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


	sysPar->emis = (double*) malloc (nr*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->emis,status,sysPar);
	sysPar->del_emit = (double*) malloc (nr*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->del_emit,status,sysPar);
	sysPar->del_inc = (double*) malloc (nr*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->del_inc,status,sysPar);

	sysPar->gstar = (double*) malloc (ng*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->gstar,status,sysPar);


	// we already set the values as they are fixed
	int ii;
	int jj;
	for (ii=0; ii<ng;ii++){
		sysPar->gstar[ii] = GFAC_H + (1.0-2*GFAC_H)/(ng-1)*( (float) (ii) );
	}

	sysPar->d_gstar = (double*) malloc (ng*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->gstar,status,sysPar);
	for (ii=0; ii<ng;ii++){
	     if ((ii==0)||(ii==(ng-1))) {
	        sysPar->d_gstar[ii] = 0.5*(sysPar->gstar[1]-sysPar->gstar[0])+GFAC_H;
	     } else {
	    	 sysPar->d_gstar[ii] = sysPar->gstar[1]-sysPar->gstar[0];
	     }
	}



	sysPar->trff = (double***) malloc(nr*sizeof(double**));
	CHECK_MALLOC_RET_STATUS(sysPar->trff,status,sysPar);
	sysPar->cosne = (double***) malloc(nr*sizeof(double**));
	CHECK_MALLOC_RET_STATUS(sysPar->cosne,status,sysPar);
	for (ii=0; ii < nr; ii++){
		sysPar->trff[ii] = (double**) malloc(ng*sizeof(double*));
		CHECK_MALLOC_RET_STATUS(sysPar->trff[ii],status,sysPar);
		sysPar->cosne[ii] = (double**) malloc(ng*sizeof(double*));
		CHECK_MALLOC_RET_STATUS(sysPar->cosne[ii],status,sysPar);
		for (jj=0; jj < ng; jj++){
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
		free(sysPar->d_gstar);

		free(sysPar->emis);
		free(sysPar->del_emit);
		free(sysPar->del_inc);

		if(sysPar->trff != NULL){
			int ii;
			for (ii=0; ii<sysPar->nr; ii++){
				if(sysPar->trff[ii] != NULL){
					int jj;
					for (jj=0; jj<sysPar->ng; jj++){
						free(sysPar->trff[ii][jj]);
					}
					free(sysPar->trff[ii]);
				}
			}
			free(sysPar->trff);
		}

		if(sysPar->cosne != NULL){
			int ii;
			for (ii=0; ii<sysPar->nr; ii++){
				if(sysPar->cosne[ii] != NULL){
					int jj;
					for (jj=0; jj<sysPar->ng; jj++){
						free(sysPar->cosne[ii][jj]);
					}
					free(sysPar->cosne[ii]);
				}
			}
			free(sysPar->cosne);
		}
		free(sysPar);
	}
}




/*** struct timeval start, end;
	long mtime, seconds, useconds;
	gettimeofday(&start, NULL);
gettimeofday(&end, NULL);

    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    mtime = ((seconds) * 1000*1000 + useconds) + 0.5;

    printf("Elapsed time: %ld micro seconds\n", mtime); ***/
