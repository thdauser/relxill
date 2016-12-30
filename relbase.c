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
rel_spec* cached_rel_spec=NULL;
double cached_int_romb_rad = -1.0;
const double cache_limit = 1e-8;


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
static void interpol_relTable(relSysPar** sysPar_inp,double a, double mu0, double rin, double rout,
		 int* status){

	// load tables
	if (relline_table==NULL){
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

	/** TODO: check if we need R_DIFF
	double r_diff = rms - ((1.0-ifac_a)*kerr_rms(tab->a[ind_a])
	       + (ifac_a)*kerr_rms(tab->a[ind_a+1]) ); **/

//	printf(" rms:%.3f, rdiff:%.3e -- %.3e\n",rms,r_diff,tab->arr[ind_a][ind_mu0]->r[tab->n_r-1]);
//	printf(" ind_a:%i ind_mu0:%i \n",ind_a,ind_mu0);

	/** get the radial grid (the radial grid only changes with A by the table definition) **/
	assert(fabs(tab->arr[ind_a][ind_mu0]->r[tab->n_r-1]
			    - tab->arr[ind_a][ind_mu0]->r[tab->n_r-1]) < 1e-6);
	for (int ii=0; ii < tab->n_r; ii++){
		cached_tab_sysPar->re[ii] = interp_lin_1d(ifac_a,
				tab->arr[ind_a][ind_mu0]->r[ii],tab->arr[ind_a+1][ind_mu0]->r[ii]);
	}

	// get the extent of the disk (indices are defined such that tab->r[ind] <= r < tab->r[ind+1]
	int ind_rmin = inv_binary_search(cached_tab_sysPar->re,tab->n_r,rin);
	int ind_rmax = inv_binary_search(cached_tab_sysPar->re,tab->n_r,rout);

	for (int ii=0; ii < tab->n_r; ii++){
		// TODO: SHOULD WE ONLY INTERPOLATE ONLY THE VALUES WE NEED??? //
		// only interpolate values where we need them
		if (ii>=ind_rmin || ii<=ind_rmax+1){
			interpol_a_mu0(ii, ifac_a, ifac_mu0, ind_a, ind_mu0, cached_tab_sysPar,tab);
		} else {  // set everything we won't need to 0 (just to be sure)
			cached_tab_sysPar->gmin[ii]=0.0;
			cached_tab_sysPar->gmax[ii]=0.0;
		    for (int jj=0; jj<tab->n_g;jj++){
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

	relSysPar* sysPar = (*sysPar_inp);
	// only need to initialize and allocat memory if not already loaded
	if (sysPar == NULL){
		sysPar = new_relSysPar(N_FRAD,tab->n_g,status);
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
	}



	return cached_relSysPar;
}

/** get new structure to store the relline spectrum (possibly for several zones)
    important note: ener has n_ener+1 number of bins **/
static rel_spec* new_rel_spec(int nzones, const int n_ener, int*status){
	rel_spec* spec = (rel_spec*) malloc(sizeof(rel_spec));
	CHECK_MALLOC_RET_STATUS(spec,status,NULL);


	spec->n_zones = nzones;
	spec->n_ener  = n_ener;

	spec->flux = (double**) malloc (spec->n_zones * sizeof(double*) );
	CHECK_MALLOC_RET_STATUS(spec->flux,status,spec);

	int ii;
	for (ii=0; ii<spec->n_ener; ii++){
		spec->flux[ii] = (double*) malloc ( n_ener * sizeof(double) );
		CHECK_MALLOC_RET_STATUS(spec,status,spec);
	}

	return spec;
}


/** get the number of zones on which we calculate the relline-spectrum **/
static int get_num_zones(int model_type){

	// set the number of zones in radial direction (1 for relline/conv model, N_ZONES for xill models)
	if (is_xill_model(model_type)){
		return N_ZONES;
	} else {
		return 1;
	}

}


/** initialize the rel_spec structure **/
static void init_rel_spec(rel_spec** spec, relParam* param, const double* ener, const int n_ener, int* status ){

	int nzones = get_num_zones(param->model_type);
	if ((*spec)==NULL){
		(*spec) = new_rel_spec(nzones,n_ener,status);
	}
	double* r_grid = get_rzone_grid(param->rin, param->rout, nzones, status);
	CHECK_STATUS_VOID(*status);
	(*spec)->ener = &ener;
	(*spec)->rgrid = r_grid;

	return;
}

static void zero_rel_spec_flux(rel_spec* spec){
	int ii; int jj;
	for (ii=0; ii<spec->n_zones;ii++){
		for (jj=0; jj<spec->n_ener;jj++){
			spec->flux[ii][jj] = 0.0;
		}
	}
}


/** Romberg Integration Routine **/
static double RombergIntegral(a,b,str_relb_func* str, int k){
  const double prec = 0.05;
  double obtprec;
  const int itermin = 0;
  const int itermax = 0;
  const int maxiter = 6;
  double t[maxiter+1][maxiter+1];

  if (itermax>maxiter) {
	  itermax=maxiter;
  }

  // check if this value has already been calculated
  double r;
  if (str->cached_relbf) {
     r = str->cache_val_relb_func[k];
  } else {
     r = RELB_FUNC(a,k);
  }

  str->cache_val_relb_func[k] = RELB_FUNC(b,k);
  str->cache_rad_relb_fun = str->re;
  // rb(k) = RELB_FUNC(b,k);

  int ii;

  double ta = (r + str->cache_val_relb_func[k]) / 2.0;
  double niter = 0;
  double pas=b-a;
  double pasm=1.0;
  t[0][0]=ta*pas;
  while ( ( niter<itermin) && ( obtprec < prec) && ( niter <= itermax)) {
	  niter++;
	  pas=pas/2.0;
	  pasm=pasm/2.0;
	  double s=ta;
	  for (ii=0; ii<niter*niter-1; ii++) {
		  s += RELB_FUNC(a+pas*ii,k);
	  }
  t[0][niter]=s*pas;
  r=1.0;
  for (ii=0; ii<niter; ii++){
     r *= 4.0;
     int jj=niter-ii;
     t[ii][jj] = (r*t[ii-1][jj+1] - t[ii-1][jj])/(r-1.0);
  }
   obtprec = fabs(t[niter][0] - t[niter-1][0])/t[niter][0];
  }

  return t[niter][0];
}


static str_relb_func* new_str_relb_func(relSysPar* sysPar, int* status){
	str_relb_func* str = (str_relb_func*) malloc(sizeof(str_relb_func));
	CHECK_MALLOC_RET(str,status,NULL);

	str->gstar = sysPar->gstar;
	str->ng = sysPar->ng;

	return str;
}
static free_str_relb_func(str_relb_func* str){
	free(str);
}

/** relat. function which we want to integrate **/
static double relb_func(double e, int k, str_relb_func* str, double line_energy){

  // get the redshift from the energy
  double eg = e/line_energy;
  double egstar = (eg-str->gmin)*str->del_g;

  // find the indices in the original g-grid, but check first if they have already been calculated
  int ind;
  if (!((egstar>=str->gstar[str->save_g_ind])&&(egstar<str->gstar[str->save_g_ind+1]))){
	  str->save_g_ind  = binary_search(egstar,str->gstar,str->ng);
  }
  ind = str->save_g_ind;

  double inte = (egstar -gstar(ind))/(gstar[ind+1] - gstar[ind]);
  double inte1=1.0-inte;
  double ftrf = inte*str->trff[ind,k] + inte1*str->trff[ind+1,k];
  double fmu0 = inte*str->cosne[ind,k] + inte1*str->cosne[ind+1,k];

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

/** function which makes an approximated integration for gstar->0/1
    this is only done within gstar=[0,H] and gstar[H,1-H]
    input:   bin_lo and bin_hi
    output:  area of the bin (= luminosity = E/dt/bin) )  **/
static double int_edge(blo,bhi,h,str_relb_func* str, line_energy){


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
    norm = norm + relb_func(gstar2ener(hex,str->gmin,str->gmax,line_energy),k);
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
			flu += RombergIntegral(lo,hi,k);
		}
	} else {
		for (k=0;k<2;k++){
			flu += relb_func((hi+lo)/2.0,k)*(hi-lo);
		}
	}

	return flu;
}

/** integrate the flux bin (see Dauser+2010, MNRAS for details) **/
static double integ_relline_bin(str_relb_func* str, double rlo0, double rhi0, double line_ener){

	double flu=0.0;

	int gblo = (rlo0/line_ener - str->gmin)*str->del_g;
	if (gblo<0.0) {
		gblo=0.0;
	} else if (gblo>1.0){
		gblo = 1.0;
	}

	int gbhi = (rhi0/line_ener - str->gmin)*str->del_g;
	if (gbhi<0.0) {
		gbhi=0.0;
	} else if (gbhi>1.0){
		gbhi = 1.0;
	}
	if (gbhi==0){
		return 0.0;
	}

	int rlo = rlo0;
	int rhi = rhi0;

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

	int calc=0;
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
	str->emis;

	str->trff = trff;
	str->cosne = cosne;

	str->cache_bin_ener = -1.0;
	str->cached_relbf = 0;
}

/** calculate the relline profile(s) for all given zones **/
str_relb_func* cached_str_relb_func = NULL;
void relline_profile(rel_spec* spec, relSysPar* sysPar, double line_ener, int* status){

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
	           spec->flux[izone,jj] +=
	        	integ_relline_bin(cached_str_relb_func,spec->ener[jj],spec->ener[jj+1],line_ener)*weight;
	     	 }
	     }
	}

}

/* the relbase function calculating the basic relativistic line shape for a given parameter setup
 * input: ener(n_ener), param
 * output: photar(n_ener)     */
void relbase(const double* ener, const int n_ener, double* photar, relParam* param, int* status){

	// check caching here?? -> free values we need to re-alloc!

	// initialize parameter values
	relSysPar* sysPar = get_system_parameters(param,status);
	CHECK_STATUS_VOID(*status);

	// get emissivity profile
	calc_emis_profile(param, sysPar, status);
	CHECK_STATUS_VOID(*status);

	// init the spectra where we store the flux
	init_rel_spec(&cached_rel_spec, param, ener, n_ener, status);
	CHECK_STATUS_VOID(*status);

	// calculate line profile
	relline_profile(cached_rel_spec, sysPar, param->lineE, status);


	// cache everything once we've calculated all the stuff
	// cached_params (!!!)
}

static free_rel_spec(rel_spec* spec){
	if (spec!=NULL){
		free(spec->ener);
		if (spec->flux!=NULL){
			int ii;
			for (ii=0; ii<spec; ii++){
				free(spec->flux[ii]);
			}
		}
		free(spec->flux);
		free(spec);
	}
}

void free_cached_tables(void){
	free_relTable(relline_table);
	free_relSysPar(cached_relSysPar);
	free_relSysPar(cached_tab_sysPar);

	free(cached_rel_param);

	free_rel_spec(cached_rel_spec);

	free_str_relb_func(cached_str_relb_func);
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

	sysPar->gstar = (double*) malloc (ng*sizeof(double));
	CHECK_MALLOC_RET_STATUS(sysPar->gstar,status,sysPar);


	// we already set the values as they are fixed
	for (int ii=0; ii<ng;ii++){
		sysPar->gstar[ii] = GFAC_H + (1.0-2*GFAC_H)/(ng-1)*( (float) (ii) );
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

		free(sysPar->emis);

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

