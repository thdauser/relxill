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

void init_par_relline(relParam* param, const double* inp_par, const int n_parameter, int* status){

	assert(n_parameter == NUM_PARAM_RELLINE);

	param->lineE = inp_par[0];
	param->emis1 = inp_par[1];
	param->emis2 = inp_par[2];
	param->rbr   = inp_par[3];
	param->a     = inp_par[4];
	param->incl  = inp_par[5];
	param->rin   = inp_par[6];
	param->rout  = inp_par[7];
	param->z     = inp_par[8];
}

void relline(const double* ener, const int n_ener, double* photar, const double* parameter, const int n_parameter, int* status){

	// fill in parameters
	relParam* param_struct = new_relParam(n_parameter,status);
	CHECK_STATUS_VOID(*status);

	init_par_relline(param_struct,parameter,n_parameter,status);
	CHECK_STATUS_VOID(*status);

	// call the function which calculates the line
	relbase(ener, n_ener, photar, param_struct,status);
	CHECK_STATUS_VOID(*status);

	free_relParam(param_struct);
}

/* get a new relbase parameter structure and initialize it */
relParam* new_relParam(int model_type, int* status){
	relParam* param = (relParam*) malloc(sizeof(relParam));
	if (param==NULL){
		RELXILL_ERROR("memory allocation failed",status);
		return NULL;
	}
	param->model_type = model_type;

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
