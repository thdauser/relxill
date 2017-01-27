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

#include "relmodels.h"


static void check_negative_radii(double* r, double a){
	if (*r<0){
		*r = -1.0*(*r)*kerr_rms(a);
	}
}

int warned_rms = 0;
static void check_parameter_bounds(relParam* param, int* status){

	// first set the Radii to positive value
	check_negative_radii(&(param->rin), param->a);
	check_negative_radii(&(param->rout), param->a);
	check_negative_radii(&(param->rbr), param->a);

	if (param->rout<=param->rin){
		printf(" *** error : Rin >= Rout not possible, please set the parameters  \n");
		*status=EXIT_FAILURE;
	}

	double rms = kerr_rms(param->a);
	if (param->rin < rms){
		if (!warned_rms){
			printf(" *** warning : Rin < ISCO, resetting Rin=ISCO; please set your limits properly \n");
		}
		param->rin = rms;
	}


	// TODO: also check Rbreak here

}


xillParam* init_par_xillver(const double* inp_par, const int n_parameter, int* status){

	// fill in parameters
	xillParam* param = new_xillParam(MOD_TYPE_XILLVER,status);
	CHECK_STATUS_RET(*status,NULL);

	assert(n_parameter == NUM_PARAM_XILLVER);

	param->gam   = inp_par[0];
	param->afe   = inp_par[1];
	param->lxi   = inp_par[2];
	param->ect   = inp_par[3];
	param->incl  = inp_par[4]; // is given in degrees !!
	param->z     = inp_par[5];

	// TODO: check parameter bounds here as well
/*	check_parameter_bounds_xillver(param,status);
	CHECK_STATUS_RET(*status,NULL); */

	return param;
}


void init_par_relxill(relParam** rel_param, xillParam** xill_param, const double* inp_par, const int n_parameter, int* status){

	// fill in parameters
	relParam* param = new_relParam(MOD_TYPE_RELXILL,EMIS_TYPE_BKN,status);
	CHECK_STATUS_VOID(*status);

	xillParam* xparam = new_xillParam(MOD_TYPE_RELXILL,status);
	CHECK_STATUS_VOID(*status);

	assert(n_parameter == NUM_PARAM_RELXILL);

	param->emis1 = inp_par[0];
	param->emis2 = inp_par[1];
	param->rbr   = inp_par[2];
	param->a     = inp_par[3];
	param->incl  = inp_par[4]*M_PI/180;
	param->rin   = inp_par[5];
	param->rout  = inp_par[6];
	param->z     = inp_par[11];

	xparam->gam   = inp_par[7];
	xparam->afe   = inp_par[8];
	xparam->lxi   = inp_par[9];
	xparam->ect   = inp_par[10];
	xparam->z     = inp_par[11];


	check_parameter_bounds(param,status);
	CHECK_STATUS_VOID(*status);

	*rel_param  = param;
	*xill_param = xparam;

	return;
}


relParam* init_par_relline(const double* inp_par, const int n_parameter, int* status){

	// fill in parameters
	relParam* param = new_relParam(MOD_TYPE_RELLINE,EMIS_TYPE_BKN,status);
	CHECK_STATUS_RET(*status,NULL);

	assert(n_parameter == NUM_PARAM_RELLINE);

	param->lineE = inp_par[0];
	param->emis1 = inp_par[1];
	param->emis2 = inp_par[2];
	param->rbr   = inp_par[3];
	param->a     = inp_par[4];
	param->incl  = inp_par[5]*M_PI/180;
	param->rin   = inp_par[6];
	param->rout  = inp_par[7];
	param->z     = inp_par[8];

	check_parameter_bounds(param,status);
	CHECK_STATUS_RET(*status,NULL);

	return param;
}

relParam* init_par_relline_lp(const double* inp_par, const int n_parameter, int* status){

	// fill in parameters
	relParam* param = new_relParam(MOD_TYPE_RELLINELP,EMIS_TYPE_LP,status);
	CHECK_STATUS_RET(*status,NULL);

	assert(n_parameter == NUM_PARAM_RELLINELP);

	param->lineE  = inp_par[0];
	param->height = inp_par[1];
	param->a      = inp_par[2];
	param->incl   = inp_par[3]*M_PI/180;
	param->rin    = inp_par[4];
	param->rout   = inp_par[5];
	param->z      = inp_par[6];
	param->gamma  = inp_par[7];

	check_parameter_bounds(param,status);
	CHECK_STATUS_RET(*status,NULL);

	return param;
}


/** shift the spectrum such that we can calculate the line for 1 keV **/
static double* shift_energ_spec_1keV(const double* ener, const int n_ener, double line_energ, double z,int* status){

	double* ener1keV = (double*) malloc((n_ener+1)*sizeof(double));
	CHECK_MALLOC_RET_STATUS(ener1keV,status,NULL);

	int ii;
	for (ii=0; ii<=n_ener; ii++){
		// ener1keV[ii] = ener[ii];
		ener1keV[ii] = ener[ii]*(z + line_energ); //TODO: need to test this
	}
	return ener1keV;
}

/** XSPEC XILLVER MODEL FUNCTION **/
void relxill(const double* ener0, const int n_ener0, double* photar, const double* parameter, const int n_parameter, int* status){

	xillParam* xill_param = NULL;
	relParam* rel_param = NULL;

	init_par_relxill(&rel_param,&xill_param,parameter,n_parameter,status);
	CHECK_STATUS_VOID(*status);

	double * ener = (double*) ener0;
	int n_ener = (int) n_ener0;
	relxill_kernel(ener, photar, n_ener, xill_param, rel_param,status);
	CHECK_STATUS_VOID(*status);

	// todo: shift spectrum accordingly (or already in xillver???)

	// test output
	save_xillver_spectrum(ener,photar,n_ener);

	free_xillParam(xill_param);
	free_relParam(rel_param);

}


/** XSPEC XILLVER MODEL FUNCTION **/
void xillver(const double* ener0, const int n_ener0, double* photar, const double* parameter, const int n_parameter, int* status){

	xillParam* param_struct = init_par_xillver(parameter,n_parameter,status);
	CHECK_STATUS_VOID(*status);

	// shift the spectrum such that we can calculate the line for 1 keV
	 /**double* ener1keV = shift_energ_spec_1keV(ener, n_ener, 0.0 , param_struct->z,status);
	 CHECK_STATUS_VOID(*status); **/

	// call the function which calculates the xillver spectrum
	xill_spec* spec = get_xillver_spectra(param_struct,status);
	CHECK_STATUS_VOID(*status);

	// rebin it to the given grid

	// =4= rebin to the input grid
	assert(spec->n_incl==1); // make sure there is only one spectrum given (for the chosen inclination)

	double * ener = (double*) ener0;
	int n_ener = (int) n_ener0;
	rebin_spectrum( ener, photar,n_ener,spec->ener,spec->flu[0],spec->n_ener);

	// test output
	save_xillver_spectrum(ener,photar,n_ener);

	free_xillParam(param_struct);

	free_xill_spec(spec);
}



/** XSPEC RELLINE MODEL FUNCTION **/
void relline(const double* ener, const int n_ener, double* photar, const double* parameter, const int n_parameter, int* status){

	relParam* param_struct = init_par_relline(parameter,n_parameter,status);
	CHECK_STATUS_VOID(*status);

	// shift the spectrum such that we can calculate the line for 1 keV
	 double* ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z,status);
	 CHECK_STATUS_VOID(*status);

	// call the function which calculates the line (assumes a line at 1keV!)
	rel_spec* spec = relbase(ener1keV, n_ener, photar, param_struct,status);
	CHECK_STATUS_VOID(*status);

	save_relline_profile(spec);

	free_relParam(param_struct);
}

/** XSPEC RELLINELP MODEL FUNCTION **/
void rellinelp(const double* ener, const int n_ener, double* photar, const double* parameter, const int n_parameter, int* status){

	relParam* param_struct = init_par_relline_lp(parameter,n_parameter,status);
	CHECK_STATUS_VOID(*status);

	// shift the spectrum such that we can calculate the line for 1 keV
	 double* ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z,status);
	 CHECK_STATUS_VOID(*status);

	// call the function which calculates the line (assumes a line at 1keV!)
	relbase(ener1keV, n_ener, photar, param_struct,status);
	CHECK_STATUS_VOID(*status);

	free_relParam(param_struct);
}


/* get a new relbase parameter structure and initialize it */
relParam* new_relParam(int model_type, int emis_type, int* status){
	relParam* param = (relParam*) malloc(sizeof(relParam));
	if (param==NULL){
		RELXILL_ERROR("memory allocation failed",status);
		return NULL;
	}
	param->model_type = model_type;
	param->emis_type  = emis_type;

	param->a = PARAM_DEFAULT;
	param->incl = PARAM_DEFAULT;
	param->emis1 = PARAM_DEFAULT;
	param->emis2 = PARAM_DEFAULT;
	param->rbr = PARAM_DEFAULT;
	param->rin = PARAM_DEFAULT;
	param->rout = PARAM_DEFAULT;
	param->lineE = PARAM_DEFAULT;
	param->z = PARAM_DEFAULT;
	param->height = PARAM_DEFAULT;
	param->gamma = PARAM_DEFAULT;

	return param;
}


/* free relbase parameter */
void free_relParam(relParam* param){
	free(param);
}



/* get a new relbase parameter structure and initialize it */
xillParam* new_xillParam(int model_type, int* status){
	xillParam* param = (xillParam*) malloc(sizeof(xillParam));
	if (param==NULL){
		RELXILL_ERROR("memory allocation failed",status);
		return NULL;
	}
	param->model_type = model_type;

	param->gam = PARAM_DEFAULT;
	param->afe = PARAM_DEFAULT;
	param->lxi = PARAM_DEFAULT;
	param->ect = PARAM_DEFAULT;
	param->incl = PARAM_DEFAULT;
	param->z = PARAM_DEFAULT;

	return param;
}


/* free relbase parameter */
void free_xillParam(xillParam* param){
	free(param);
}
