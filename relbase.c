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


// new CACHE routines
cnode* cache_relbase = NULL;
cnode* cache_syspar  = NULL;


/** global parameters, which can be used for several calls of the model */
relTable* relline_table=NULL;
relSysPar* cached_tab_sysPar=NULL;

/** caching parameters **/
relParam* cached_rel_param=NULL;
xillParam* cached_xill_param=NULL;


int save_1eV_pos = 0;
double cached_int_romb_rad = -1.0;
//const double cache_limit = 1e-8;

const double ener_xill_norm_lo = 0.1;
const double ener_xill_norm_hi = 1000;
const int n_ener_xill = 3000;
double* ener_xill = NULL;

const int n_ener_std = N_ENER_CONV;
double* ener_std = NULL;

int c_num_zones=0;

specCache* spec_cache = NULL;

// precision to calculate gstar from [H:1-H] instead of [0:1]
const double GFAC_H = 5e-3;

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

	double r1=1.0/sqrt(rout);
	double r2=1.0/sqrt(rin);
	int ii;
	for (ii=0; ii<sysPar->nr; ii++){
		sysPar->re[ii] = ((double) (ii) )*(r2-r1)/(sysPar->nr-1)+r1;
		sysPar->re[ii] = pow(1.0/(sysPar->re[ii]),2);
		assert(sysPar->re[ii]>1.0);
	}
	return;
}

/* function interpolating the rel table values for rin,rout,mu0,incl   */
static relSysPar* interpol_relTable(double a, double mu0, double rin, double rout,
		 int* status){

	// load tables
	if (relline_table==NULL){
		print_version_number(status);
		CHECK_STATUS_RET(*status,NULL);
		read_relline_table(RELTABLE_FILENAME,&relline_table,status);
		CHECK_STATUS_RET(*status,NULL);
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
		CHECK_STATUS_RET(*status,NULL);
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

	//  need to initialize and allocate memory
	relSysPar* sysPar = new_relSysPar(N_FRAD,tab->n_g,status);
	CHECK_STATUS_RET(*status,NULL);
	get_fine_radial_grid(rin,rout,sysPar);

	/** we do not have rmax=1000.0 in the table, but just values close to it so let's do this trick**/
	double rmax = 1000.0;
	if (cached_tab_sysPar->re[ind_rmax] < rmax && cached_tab_sysPar->re[ind_rmax]*1.01 > rmax ){
		cached_tab_sysPar->re[ind_rmax] = rmax;
	}

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
					CHECK_STATUS_RET(*status,NULL);
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
			CHECK_STATUS_RET(*status,NULL);
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

	return sysPar;
}


/**  calculate all relativistic system parameters, including interpolation
 *   of the rel-table, and the emissivity; caching is implemented
 *   Input: relParam* param   Output: relSysPar* system_parameter_struct
 */


/* function to get the system parameters */
static relSysPar* calculate_system_parameters(relParam* param, int* status){

	CHECK_STATUS_RET(*status,NULL);

	// only re-do the interpolation if rmin,rmax,a,mu0 changed
	// or if the cached parameters are NULL

	double mu0 = cos(param->incl);
	relSysPar* sysPar = interpol_relTable(param->a,mu0,param->rin,param->rout,status);
	CHECK_STATUS_RET(*status,NULL);


	if (param->limb!=0){
		sysPar->limb_law = param->limb;
	}


	// if we are in LP Geometry, we also need primary source emission angle to hit Rin and Rout
	if (param->emis_type==EMIS_TYPE_LP){
		get_ad_del_lim(param, sysPar, status);
	}


	// get emissivity profile
	calc_emis_profile(sysPar->emis, sysPar->del_emit, sysPar->del_inc, sysPar->re, sysPar->nr, param, status);

	if (*status!=EXIT_SUCCESS){
		RELXILL_ERROR("failed to calculate the system parameters",status);
	}

	return sysPar;
}


relSysPar* get_system_parameters(relParam* param, int* status){

	CHECK_STATUS_RET(*status, NULL);

	inpar* sysinp = set_input_syspar(param,status);
	CHECK_STATUS_RET(*status,NULL);

	cache_info* ca_info = cli_check_cache(cache_syspar, sysinp, check_cache_syspar, status);
	CHECK_STATUS_RET(*status,NULL);

	relSysPar* sysPar = NULL;
	if (ca_info->syscache==1){
		// system parameter values are cached, so we can take it from there
		sysPar = ca_info->store->data->relSysPar;
		if (is_debug_run()){
			printf(" DEBUG:  SYSPAR-Cache: re-using calculated values\n");
		}
	} else {
		// NOT CACHED, so we need to calculate the system parameters
		sysPar = calculate_system_parameters(param, status);

		// now add (i.e., prepend) the current calculation to the cache
		set_cache_syspar(&cache_syspar, param, sysPar,status);

		if (is_debug_run() && *status==EXIT_SUCCESS){
			printf(" DEBUG:  The count of the SYSPAR-Cache  is  %i \n",cli_count_elements(cache_syspar));
		}
	}

	free(ca_info);
	free(sysinp);

	CHECK_RELXILL_DEFAULT_ERROR(status);

	// make a sanity check for now
	if (*status == EXIT_SUCCESS){
		assert(cache_syspar!=NULL);
		assert(sysPar!=NULL);
	}

	return  sysPar;
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

	CHECK_STATUS_VOID(*status);

	/** in case of the relxill-LP model multiple zones are used **/
	int nzones = param->num_zones; ///  get_num_zones(param->model_type, param->emis_type, ION_GRAD_TYPE_CONST);

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

	get_rzone_grid(param->rin, param->rout, (*spec)->rgrid, nzones, param->height, status);

	CHECK_RELXILL_DEFAULT_ERROR(status);

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
	str->limb_law = 0;

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

  double val = pow(eg,3)/((str->gmax-str->gmin)*sqrt(egstar - egstar*egstar))*ftrf*str->emis;

  /** isotropic limb law by default (see Svoboda (2009)) **/
  if (str->limb_law==0){
	  return val;
  } else {
	  double fmu0 = inte*str->cosne[ind][k] + inte1*str->cosne[ind+1][k];
	  double limb = 1.0;
	  if (str->limb_law==1) { //   !Laor(1991)
		  limb = (1.0 + 2.06*fmu0);
	  }else if (str->limb_law == 2) {  //  !Haardt (1993)
		  limb = log(1.0+1.0/fmu0);
	  } else {
		  return val;
	  }
	  return val*limb;
  }
}

/** Romberg Integration Routine **/
static double romberg_integration(double a,double b,int k, str_relb_func* str, double line_ener){
  const double prec = 0.02;
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


static void free_str_relb_func(str_relb_func** str){
	if (*str!=NULL){
		free(*str);
		*str=NULL;
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

	/** We are doing a trick here: for the "red wing" of the line, one can show that a simple trapez integration is enough,
	 * as the profile is very smooth. In order to avoid problems for narrow, double peaked lines, the red wing is defined
	 * to start at 0.95*line_energy **/
	int k;
	if (lo>=line_energy*0.95){
		for (k=0;k<2;k++){
			flu += romberg_integration(lo,hi,k,str,line_energy);
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
		if ( (fabs(rlo - str->cache_bin_ener) < CACHE_LIMIT) && (fabs(str->re - str->cache_rad_relb_fun) < CACHE_LIMIT)) {
			str->cached_relbf = 1;
		} else {
			str->cached_relbf = 0;
		}
		flu = flu  + int_romb(rlo,rhi,str,line_ener);
	}

	return flu;
}

static void	set_str_relbf(str_relb_func* str, double re, double gmin, double gmax, double** trff,
		double** cosne, double emis, int limb_law){
	str->re = re;
	str->gmin = gmin;
	str->gmax = gmax;
	str->del_g = 1./(gmax-gmin);
	str->emis = emis;

	str->trff = trff;
	str->cosne = cosne;

	str->cache_bin_ener = -1.0;
	str->cache_rad_relb_fun = -1.0;
	str->cached_relbf = 0;

	str->limb_law = limb_law;

	str->save_g_ind = 0;
}


/** function to properly re-normalize the relline_profile **/
static void renorm_relline_profile(rel_spec* spec, relParam* rel_param, int* status){

	CHECK_STATUS_VOID(*status);

	// normalize to 'cts/bin'
	int ii; int jj;
	double sum = 0.0;
	for (ii=0; ii<spec->n_zones; ii++){
		for (jj=0; jj<spec->n_ener; jj++){
			spec->flux[ii][jj] /= 0.5*( spec->ener[jj] + spec->ener[jj+1] );
			sum += spec->flux[ii][jj];
		}
	}

	/** only renormalize if not the relxill model or not a lamp post model **/

	if (do_renorm_model(rel_param)) {
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

	CHECK_STATUS_VOID(*status);

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

	       	// set the current parameters in a cached structure (and reset some values) [optimizes speed]
	       	set_str_relbf(cached_str_relb_func,
	       			sysPar->re[ii], sysPar->gmin[ii],sysPar->gmax[ii],
	       			sysPar->trff[ii],sysPar->cosne[ii],
					sysPar->emis[ii], sysPar->limb_law );

	        /** INTEGRATION
	         *   [remember: defintion of Xillver/Relxill is 1/2 * Speith Code]
	         *   [remember: trapez integration returns just r*dr*PI (full integral is over dA/2)]
	         *   [ -> in the end it's weigth=PI*r*dr/2 ]
	         */
	       	double weight = trapez_integ_single(sysPar->re,ii,sysPar->nr)/2;

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

	/** we need to free the structure as it points to the currently cached sysPar structure
	 	which is freed if the cache is full and therefore causes "invalid reads" **/
	free_str_relb_func(&cached_str_relb_func);

	CHECK_RELXILL_DEFAULT_ERROR(status);

}

/** convolve the (bin-integrated) spectra f1 and f2 (which need to have a certain binning)
 *  fout: gives the output
 *  f1 input (reflection) specrum
 *  f2 filter
 *  ener has length n+1 and is the energy array
 * **/
void fft_conv_spectrum(double* ener, double* fxill, double* frel, double* fout, int n,
		int re_rel, int re_xill, int izone, int* status){

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
	double xcomb[n]; double ycomb[n];

	/**********************************************************************/
	/** cache either the relat. or the xillver part, as only one of the
	 * two changes most of the time (reduce time by 1/3 for convolution) **/
	/**********************************************************************/

	/** #1: for the xillver part **/
	if (re_xill){
		for (ii=0; ii<n; ii++){
			spec_cache->fft_xill[izone][0][ii] = fxill[ii]/(ener[ii+1]-ener[ii]);
			spec_cache->fft_xill[izone][1][ii] = 0.0;
		}
		FFT_R2CT(1, m, spec_cache->fft_xill[izone][0], spec_cache->fft_xill[izone][1]);
	}
	double* x1 = spec_cache->fft_xill[izone][0];
	double* y1 = spec_cache->fft_xill[izone][1];

	/** #2: for the relat. part **/
 	if (re_rel){
		for (ii=0; ii<n; ii++){
			irot = (ii-save_1eV_pos+n) % n ;
			spec_cache->fft_rel[izone][0][irot ] = frel[ii]/(ener[ii+1]-ener[ii]);
			spec_cache->fft_rel[izone][1][ii] = 0.0;
		}
		FFT_R2CT(1, m, spec_cache->fft_rel[izone][0], spec_cache->fft_rel[izone][1]);
	}
	double* x2 = spec_cache->fft_rel[izone][0];
	double* y2 = spec_cache->fft_rel[izone][1];


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
	}

}


static void print_angle_dist(rel_cosne* spec, int izone){


	FILE* fp =  fopen ( "test_angle_dist.dat","w+" );
	int ii;
	for(ii=0; ii<spec->n_cosne; ii++ ){
		fprintf(fp," %e %e \n",spec->cosne[ii],
				spec->dist[izone][ii]);
	}
	if (fclose(fp)) exit(1);

}


static void calc_xillver_angdep(double* xill_flux, xill_spec* xill_spec,
		double* cosne, double* dist, int n_incl, int* status){

	CHECK_STATUS_VOID(*status);

	int ii; int jj;
	for (ii=0; ii<xill_spec->n_ener;ii++){
		xill_flux[ii] = 0.0;
	}

	for (ii=0; ii<xill_spec->n_incl;ii++){
		for (jj=0; jj<xill_spec->n_ener;jj++){
			xill_flux[jj] += dist[ii]*xill_spec->flu[ii][jj];
		}
	}

}

static double get_rzone_energ_shift(relParam* param, double rad, double del_emit){

	// del_emit: any value (for beta!=0 is okay, as it does not matter for the energy shift!)
	return gi_potential_lp(rad,param->a,param->height,param->beta,del_emit);
}

/** check if the complete spectrum is cached and if the energy grid did not change
 *  -> additionally we adapte the energy grid to the new dimensions and energy values
 *     if it's not the case yet
 */
static int is_all_cached(specCache* spec_cache,int n_ener_inp,double* ener_inp, int recompute_xill, int recompute_rel, int* status){

	if ( (recompute_xill+recompute_rel) == 0 && ( spec_cache->out_spec != NULL )) {
		/** first need to check if the energy grid did change (might happen) **/
		if (spec_cache->out_spec->n_ener != n_ener_inp){
			return 0;
		}
		int ii;
		// we know that n_ener
		for (ii=0; ii<n_ener_inp; ii++){
			if ( fabs(spec_cache->out_spec->ener[ii] - ener_inp[ii]) > 1e-4){
				int jj;
				for (jj=0; jj<n_ener_inp; jj++){
					spec_cache->out_spec->ener[jj] = ener_inp[jj];
				}
				return 0;
			}
		}
		return 1;
	} else {
		return 0;
	}

}


static void check_caching_relxill(relParam* rel_param, xillParam* xill_param, int* re_rel, int* re_xill){


	/** always re-compute if the number of zones changed **/
	if (cached_rel_param != NULL){

		if (rel_param->num_zones != cached_rel_param->num_zones){
			if ( is_debug_run() ){
				printf("  *** warning :  the number of radial zones was changed from %i to %i \n",
						rel_param->num_zones, cached_rel_param->num_zones);
			}
			*re_rel  = 1;
			*re_xill = 1;
			return;
		}
	} else {
		*re_rel  = 1;
		*re_xill = 1;
		return;
	}

	/** did any of the relat. parameters change?  **/
	*re_rel  = redo_relbase_calc(rel_param, cached_rel_param);

	*re_xill = redo_xillver_calc(rel_param, xill_param, cached_rel_param, cached_xill_param);

	return;
}

/** renorm a model (flu) to have the same flux as another model (flu)
 *  (bin-integrated flux, same energy grid!) **/
static void renorm_model(double* flu0, double* flu, int nbins){

	double sum_inp = 0.0;
	double sum_out = 0.0;
	int ii;
	for (ii=0; ii<nbins;ii++){
		sum_inp += flu0[ii];
		sum_out += flu[ii];
	}
	for (ii=0; ii<nbins;ii++){
		flu[ii] *= sum_inp / sum_out;
	}


}


static void renorm_xill_spec(double* spec , int n, double lxi, double dens){

	for (int ii=0; ii<n; ii++){
		spec[ii] /= pow(10,lxi);
		if (fabs(dens - 15) > 1e-6 ){
			spec[ii] /= pow(10,dens - 15);
		}
	}
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
	init_specCache(&spec_cache,status);
	CHECK_STATUS_VOID(*status);
	fft_conv_spectrum(ener, rebin_flux, rel_profile->flux[0], conv_out,  n_ener,
			1,1,0,status);
	CHECK_STATUS_VOID(*status);

	// need to renormalize the convolution? (not that only LP has a physical norm!!)
	if (! do_not_normalize_relline()){
		renorm_model(rebin_flux,conv_out,n_ener);
	}

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
	CHECK_STATUS_VOID(*status);

	// in case of an ionization gradient, we need to update the number of zones
	if (is_iongrad_model(rel_param->model_type, xill_param->ion_grad_type) ){
		rel_param->num_zones = get_num_zones(rel_param->model_type, rel_param->emis_type, xill_param->ion_grad_type);
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

	/** be careful, as xill_param->ect can get over-written so save the following values**/
	double ecut0 = xill_param->ect;

	/** note that in case of the nthcomp model Ecut is in the frame of the primary source
	    but for the bkn_powerlaw it is given in the observer frame */
	double ecut_primary = 0.0;
	if (xill_param->prim_type == PRIM_SPEC_ECUT ){
		ecut_primary = ecut0 * ( 1 + grav_redshift(rel_param) );
	} else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP ){
		ecut_primary = ecut0 ;
	}

	int recompute_xill=1;
	int recompute_rel=1;
	check_caching_relxill(rel_param,xill_param,&recompute_rel,&recompute_xill);

	init_specCache(&spec_cache,status);
	CHECK_STATUS_VOID(*status);

	/** is both already cached we can see if we can simply use the output flux value **/
	if ( is_all_cached(spec_cache,n_ener_inp,ener_inp,recompute_xill,recompute_rel,status) ){
		CHECK_STATUS_VOID(*status);
		for (ii=0; ii<n_ener_inp; ii++){
			spec_inp[ii] = spec_cache->out_spec->flux[ii];
		}

    /** if NOT, we need to do a whole lot of COMPUTATIONS **/
	} else {
		CHECK_STATUS_VOID(*status);

		/* *** first, stored the parameters for which we are calculating **/
		set_cached_xill_param(xill_param , &cached_xill_param, status);
		set_cached_rel_param(rel_param , &cached_rel_param, status);
		c_num_zones = rel_param->num_zones;

		/* calculate the relline profile **/

		rel_spec* rel_profile = relbase(ener, n_ener, rel_param, xill_tab, status);
		CHECK_STATUS_VOID(*status);

		/* init the xillver spectrum structure **/
		xill_spec* xill_spec = NULL;


		// currently only working for the LP version (as relxill always has just 1 zone)
		ion_grad* ion = NULL;
		if (is_iongrad_model(rel_param->model_type, xill_param->ion_grad_type) ){

			// make sure the number of zones is correctly set:
			assert(rel_param->num_zones == get_num_zones(rel_param->model_type, rel_param->emis_type, xill_param->ion_grad_type));

			ion = calc_ion_gradient(rel_param, xill_param->lxi, xill_param->ion_grad_index, xill_param->ion_grad_type,
					rel_profile->rgrid, rel_profile->n_zones, status);
			CHECK_STATUS_VOID(*status);
		}


		int jj;
		for (ii=0; ii<rel_profile->n_zones;ii++){
			assert(spec_cache!=NULL);

			// now calculate the reflection spectra for each zone (using the angular distribution)
			assert(rel_profile->rel_cosne != NULL);
			if (is_debug_run()==1){
				print_angle_dist(rel_profile->rel_cosne,ii);
			}

			// get the energy shift in order to calculate the proper Ecut value (if nzones>1)
			// (the latter part of the IF is a trick to get the same effect as NZONES=1 if during a running
			//  session the number of zones is reset)
			if (rel_profile->n_zones==1) {
//			if (rel_profile->n_zones==1 || get_num_zones(rel_param->model_type,rel_param->emis_type, xill_param->ion_grad_type)==1){
				xill_param->ect = ecut0 ;
			} else {
				// choose the (linear) middle of the radial zone
				double rzone = 0.5*(rel_profile->rgrid[ii]+rel_profile->rgrid[ii+1]);
				double del_emit = 0.0;  // only relevant if beta!=0 (doppler boosting)
				if (ion!=NULL){
					assert(ion->nbins==rel_profile->n_zones);
					del_emit = ion->del_emit[ii];
				}
				xill_param->ect = ecut_primary * get_rzone_energ_shift(rel_param,rzone, del_emit ) ;
			}

			// if we have an ionization gradient defined, we need to set the xlxi to the value of the current zone
			if (ion!=NULL){
				xill_param->lxi = ion->lxi[ii];
			}

			// call the function which calculates the xillver spectrum
			//  - always need to re-compute if we have an ionization gradient, TODO: better caching here
//			if (recompute_xill || ion!=NULL){
			if (recompute_xill){
				if (spec_cache->xill_spec[ii]!=NULL){
					free_xill_spec(spec_cache->xill_spec[ii]);
				}
				spec_cache->xill_spec[ii] = get_xillver_spectra(xill_param,status);
			}
			xill_spec = spec_cache->xill_spec[ii];

			// get angular distribution
			calc_xillver_angdep(xill_angdist_inp,xill_spec,rel_profile->rel_cosne->cosne,
					rel_profile->rel_cosne->dist[ii],rel_profile->rel_cosne->n_cosne,status);

			rebin_spectrum(ener,xill_flux,n_ener,
					xill_spec->ener, xill_angdist_inp, xill_spec->n_ener );


			/** convolve the spectrum **
			 * (important for the convolution: need to recompute fft for xillver
			 *  always if rel changes, as the angular distribution changes !!)
			 */
			fft_conv_spectrum(ener, xill_flux, rel_profile->flux[ii], conv_out,  n_ener,
					recompute_rel, 1, ii, status);
			CHECK_STATUS_VOID(*status);



			double test_sum_relline = 0.0;
			double test_sum_xillver = 0.0;
			double test_sum_relxill = 0.0;
			for (jj=0; jj<n_ener;jj++){
				if (ener[jj] > EMIN_XILLVER && ener[jj+1] < EMAX_XILLVER  ){
					test_sum_relline += rel_profile->flux[ii][jj];
					test_sum_relxill += conv_out[jj];
					test_sum_xillver += xill_flux[jj];
				}
			}

			// rebin to the output grid
			rebin_spectrum(ener_inp,single_spec_inp,n_ener_inp, ener, conv_out, n_ener);


			/** avoid problems where no relxill bin falls into an ionization bin **/
			if (test_sum_relline < 1e-12){
				continue;
			}

			// renorm the spectrum (such that it is independent of xi and density)
			renorm_xill_spec(single_spec_inp, n_ener_inp, xill_param->lxi, xill_param->dens);

			// add it to the final output spectrum
			for (jj=0; jj<n_ener_inp;jj++){
				spec_inp[jj] += single_spec_inp[jj]*test_sum_relline*test_sum_xillver/test_sum_relxill;
			}

			if (is_debug_run() && rel_profile->n_zones <= 10 ){
				char* vstr;
				double test_flu[n_ener_inp];
				for (jj=0; jj<n_ener_inp;jj++){
					test_flu[jj] = single_spec_inp[jj]*test_sum_relline*test_sum_xillver/test_sum_relxill;
				}
				if (asprintf(&vstr, "test_relxill_spec_zones_%03i.dat", ii+1) == -1){
					RELXILL_ERROR("failed to get filename",status);
				}
				save_xillver_spectrum(ener_inp,test_flu,n_ener_inp,vstr);
			}

		} /**** END OF LOOP OVER RADIAL ZONES [ii] *****/

		/** important: set the cutoff energy value back to its original value **/
		xill_param->ect = ecut0;

		/** free the ionization parameter structure **/
		free_ion_grad(ion);

		/** initialize the cached output spec array **/
		if ((spec_cache->out_spec != NULL ) ) {
			if ( spec_cache->out_spec->n_ener != n_ener_inp){
				free_out_spec(spec_cache->out_spec);
				spec_cache->out_spec = init_out_spec(n_ener_inp,ener_inp,status);
				CHECK_STATUS_VOID(*status);
			}
		} else {
			spec_cache->out_spec = init_out_spec(n_ener_inp,ener_inp,status);
			CHECK_STATUS_VOID(*status);
		}

			// lastely, we make the spectrum normalization independent of the ionization parameter
		for (ii=0; ii<n_ener_inp; ii++){
			/**		spec_inp[ii] /= pow(10,xill_param->lxi);  // RE-NORMALIZATION is now done before
			if (fabs(xill_param->dens - 15) > 1e-6 ){
				spec_inp[ii] /= pow(10,xill_param->dens - 15);
			} **/
			spec_cache->out_spec->flux[ii] = spec_inp[ii];
		}
	} /************* END OF THE HUGE COMPUTATION ********************/



	/** add a primary spectral component and normalize according to the given refl_frac parameter**/
	add_primary_component(ener_inp,n_ener_inp,spec_inp,rel_param,xill_param, status);


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

	/** note that in case of the nthcomp model Ecut is in the frame of the primary source
	    but for the bkn_powerlaw it is given in the observer frame */

	/** 1 **  calculate the primary spectrum  **/
	if (xill_param->prim_type == PRIM_SPEC_ECUT ){

		/** IMPORTANT: defintion of Ecut is ALWAYS in the frame of the observer by definition **/
		/**    (in case of the nthcomp primary continuum ect is actually kte ) **/
		double ecut_rest = xill_param->ect;

		for (ii=0; ii<n_ener_xill; ii++){
			pl_flux_xill[ii] = exp(1.0/ecut_rest) *
		             pow(0.5*(ener_xill[ii]+ener_xill[ii+1]),-xill_param->gam) *
		             exp( -0.5*(ener_xill[ii]+ener_xill[ii+1])/ecut_rest) *
		             (ener_xill[ii+1] - ener_xill[ii]);
		}

	} else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP) {
		double nthcomp_param[5];
		/** important, kTe is given in the primary source frame, so we have to add the redshift here
		 *     however, only if the REL model **/
		double z=0.0;
		if (rel_param != NULL &&rel_param->emis_type==EMIS_TYPE_LP ){
			z = grav_redshift(rel_param);
		}
		get_nthcomp_param(nthcomp_param, xill_param->gam, xill_param->ect, z);
		c_donthcomp(ener_xill, n_ener_xill, nthcomp_param, pl_flux_xill);
	} else {
		RELXILL_ERROR("trying to add a primary continuum to a model where this does not make sense (should not happen!)",status);
		return;
	}



	/** 2 **  get the normalization of the spectrum with respect to xillver **/
	/** everything is normalized to 10^15 cm^3 **/
	double norm_xill = 1e15 / 4.0 / M_PI;
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

		// should be cached, as it has been calculated before
		relSysPar* sysPar = get_system_parameters(rel_param, status);

		/** 3 **  calculate predicted reflection fraction and check if we want to use this value **/
		lpReflFrac* struct_refl_frac = calc_refl_frac(sysPar, rel_param,status);
		CHECK_STATUS_VOID(*status);

		if ((xill_param->fixReflFrac==1)||(xill_param->fixReflFrac==2)) {
			/** set the reflection fraction calculated from the height and
			 *  spin of the primary source, in this case for the physical
			 *  value from Rin to Rout          						 */
			xill_param->refl_frac = struct_refl_frac->refl_frac;
		} else {

		}

		/** 4 ** and apply it to primary and reflected spectra **/
		if (rel_param->emis_type == EMIS_TYPE_LP) {
			double g_inf = sqrt( 1.0 - ( 2*rel_param->height /
					(rel_param->height*rel_param->height + rel_param->a*rel_param->a)) );



			 // calculate the fraction of photons hitting the accretion disk
			 /** if the user sets the refl_frac parameter manually, we need to calculate the ratio
			  *  to end up with the correct normalization
			  */
			double norm_fac_refl = (fabs(xill_param->refl_frac))/struct_refl_frac->refl_frac;

 			double prim_fac = struct_refl_frac->f_inf / 0.5 * pow(g_inf,xill_param->gam);

			for (ii=0; ii<n_ener; ii++) {
				pl_flux[ii] *= norm_pl * prim_fac;
				flu[ii] *= norm_fac_refl;
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

			printf("For a = %.3f and h = %.2f rg", rel_param->a,rel_param->height);
			if (is_iongrad_model(rel_param->model_type, xill_param->ion_grad_type) || rel_param->beta>1e-6){
				printf(" and beta=%.3f v/c", rel_param->beta);
			}
			printf(": \n - reflection fraction  %.3f \n - reflection strength is: %.3f \n",	struct_refl_frac->refl_frac,sum/sum_pl);
			printf(" - %.2f%% of the photons are falling into the black hole\n", struct_refl_frac->f_bh*100);
			printf(" - gravitational redshift from the observer to the primary source is %.3f\n",grav_redshift(rel_param) );
			// TODO: add the movement of the source here
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

	if (spec==NULL) return;

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

int comp_xill_param(xillParam* cpar, xillParam* par){
	if (comp_single_param_val(par->afe,cpar->afe)) return 1;
	if (comp_single_param_val(par->dens,cpar->dens)) return 1;
	if (comp_single_param_val(par->ect,cpar->ect)) return 1;
	if (comp_single_param_val(par->gam,cpar->gam)) return 1;
	if (comp_single_param_val(par->lxi,cpar->lxi)) return 1;
	if (comp_single_param_val(par->z,cpar->z)) return 1;

	if (comp_single_param_val( (double) par->prim_type, (double) cpar->prim_type)) return 1;
	if (comp_single_param_val( (double) par->model_type, (double) cpar->model_type)) return 1;

	if (comp_single_param_val( par->ion_grad_index, cpar->ion_grad_index)) return 1;
	if (comp_single_param_val( (double) par->ion_grad_type, (double) cpar->ion_grad_type)) return 1;

	return 0;
}

/* check if values, which need a re-computation of the relline profile, have changed */
int redo_xillver_calc(relParam* rel_param, xillParam* xill_param, relParam* ca_rel_param, xillParam* ca_xill_param){

	int redo = 1;

	if ((ca_rel_param==NULL)|| (ca_xill_param==NULL)  ){
	} else {
		redo = comp_xill_param(ca_xill_param,xill_param);

		/** did spin or h change (means xillver needs to be re-computed as well, due to Ecut) **/
		if (comp_single_param_val(rel_param->a,ca_rel_param->a) ||
				comp_single_param_val(rel_param->height,ca_rel_param->height)	){
			redo = 1;
		}

	}

	return redo;
}

int redo_relbase_calc(relParam* rel_param, relParam* ca_rel_param){

	int redo = 1;
	int not_redo = 0;


	if (comp_rel_param(ca_rel_param, rel_param)){
		return redo;
	}

	return not_redo;

}

/* the relbase function calculating the basic relativistic line shape for a given parameter setup
 * (assuming a 1keV line, by a grid given in keV!)
 * input: ener(n_ener), param
 * optinal input: xillver grid
 * output: photar(n_ener)     */

rel_spec* relbase(double* ener, const int n_ener, relParam* param, xillTable* xill_tab, int* status){

	CHECK_STATUS_RET(*status,NULL);

	inpar* inp = set_input(ener,n_ener,param,NULL, status);

	// check caching here and also re-set the cached parameter values
	cache_info* ca_info = cli_check_cache(cache_relbase, inp, check_cache_relpar, status);


	// set a pointer to the spectrum
	rel_spec* spec = NULL;

	if ( is_relbase_cached(ca_info)==0 ) {


		// initialize parameter values
		relSysPar* sysPar = get_system_parameters(param,status);

		if (is_debug_run() && sysPar!=NULL){
			save_emis_profile(sysPar->re, sysPar->emis, sysPar->nr);
		}

		// init the spectra where we store the flux
		init_rel_spec(&spec, param, xill_tab, &ener, n_ener, status);

		// calculate line profile (returned units are 'cts/bin')
		relline_profile(spec, sysPar, status);

		// normalize it and calculate the angular distribution (if necessary)
		renorm_relline_profile(spec,param,status);

		if (is_debug_run()){
			save_relline_profile(spec);
		}

		// last step: store parameters and cached rel_spec (this prepends a new node to the cache)
		set_cache_relbase(&cache_relbase,param,spec, status);
		if (is_debug_run() && *status==EXIT_SUCCESS){
			printf(" DEBUG:  Adding new RELBASE eval to cache; the count is %i \n",cli_count_elements(cache_relbase));
		}
	}

	else {
		if (is_debug_run()){
			printf(" DEBUG:  RELBASE-Cache: re-using calculated values\n");
		}
		spec = ca_info->store->data->relbase_spec;
	}

	// free the input structure
	free(inp);
	free(ca_info);

	CHECK_RELXILL_DEFAULT_ERROR(status);

	return spec;
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
	free_relSysPar(cached_tab_sysPar);
	free_cached_lpTable();

	free_cached_xillTable();

	free(cached_rel_param);
	free(cached_xill_param);

	free_str_relb_func(&cached_str_relb_func);

	free_specCache();

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

	sysPar->limb_law = 0;

	sysPar->del_ad_risco=0;
	sysPar->del_ad_rmax=M_PI/2;

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


void free_fft_cache(double*** sp, int n1, int n2){

	int ii; int jj;
	if (sp !=NULL){
		for (ii=0; ii<n1; ii++){
			if (sp[ii] != NULL){
				for (jj=0; jj<n2; jj++){
					free(sp[ii][jj]);
				}
			}
			free(sp[ii]);
		}
		free(sp);
	}

}


static specCache* new_specCache(int n_cache, int n_ener, int* status){



	specCache* spec = (specCache*) malloc(sizeof(specCache));
	CHECK_MALLOC_RET_STATUS(spec,status,NULL);

	spec->n_cache = n_cache;
	spec->nzones = 0;
	spec->n_ener = n_ener;

	spec->fft_xill = (double***) malloc (sizeof(double**)*n_cache);
	CHECK_MALLOC_RET_STATUS(spec->fft_xill,status,NULL);

	spec->fft_rel = (double***) malloc (sizeof(double**)*n_cache);
	CHECK_MALLOC_RET_STATUS(spec->fft_rel,status,NULL);

	spec->xill_spec = (xill_spec**) malloc(sizeof(xill_spec*)*n_cache);
	CHECK_MALLOC_RET_STATUS(spec->xill_spec,status,NULL);


	int ii; int jj;
	int m = 2;
	for (ii=0; ii<n_cache; ii++){
		spec->fft_xill[ii] = (double**) malloc (sizeof(double*)*m);
		CHECK_MALLOC_RET_STATUS(spec->fft_xill[ii],status,NULL);
		spec->fft_rel[ii] = (double**) malloc (sizeof(double*)*m);
		CHECK_MALLOC_RET_STATUS(spec->fft_rel[ii],status,NULL);

		for (jj=0; jj<m; jj++){
			spec->fft_xill[ii][jj] = (double*) malloc (sizeof(double)*n_ener);
			CHECK_MALLOC_RET_STATUS(spec->fft_xill[ii][jj],status,NULL);
			spec->fft_rel[ii][jj] = (double*) malloc (sizeof(double)*n_ener);
			CHECK_MALLOC_RET_STATUS(spec->fft_rel[ii][jj],status,NULL);
		}

		spec->xill_spec[ii] = NULL;

	}

	spec->out_spec = NULL;

	return spec;
}


out_spec* init_out_spec(int n_ener, double* ener, int* status){

	out_spec* spec = (out_spec*) malloc( sizeof(out_spec));
	CHECK_MALLOC_RET_STATUS(spec,status,NULL);

	spec->n_ener = n_ener;
	spec->ener = (double*) malloc(sizeof(double)*n_ener) ;
	CHECK_MALLOC_RET_STATUS(spec->ener,status,NULL);
	spec->flux = (double*) malloc(sizeof(double)*n_ener) ;
	CHECK_MALLOC_RET_STATUS(spec->flux,status,NULL);

	int ii;
	for (ii=0; ii<n_ener; ii++){
		spec->ener[ii] = ener[ii];
		spec->flux[ii] = 0.0;
	}

	return spec;
}

void free_out_spec(out_spec* spec){
	if (spec!=NULL){
		free(spec->ener);
		free(spec->flux);
		free(spec);
		spec = NULL;
	}
}


void free_specCache(void){

	int ii;
	int m = 2;
	if (spec_cache != NULL){
		if (spec_cache->xill_spec != NULL){
			for (ii=0; ii<spec_cache->n_cache; ii++){
				if (spec_cache->xill_spec[ii] != NULL){
					free_xill_spec(spec_cache->xill_spec[ii]);
				}
			}
			free(spec_cache->xill_spec);
		}


		if (spec_cache->fft_xill != NULL){
			free_fft_cache(spec_cache->fft_xill,spec_cache->n_cache,m);
		}

		if (spec_cache->fft_rel != NULL){
			free_fft_cache(spec_cache->fft_rel,spec_cache->n_cache,m);
		}

		free_out_spec(spec_cache->out_spec);

	}


	free(spec_cache);

}

void init_specCache(specCache** spec, int* status){

	if ((*spec)==NULL){
		(*spec) = new_specCache(N_ZONES_MAX,N_ENER_CONV,status);
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
