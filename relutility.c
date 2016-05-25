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

#include "relutility.h"

void relxill_error(const char* const func, const char* const msg, int* status){
	*status = EXIT_FAILURE;
	printf(" *** error in relxill (%s): %s!\n", func, msg);
}

void check_relxill_error(const char* const func, const char* const msg, int* status){
	if (*status!=EXIT_SUCCESS){
		*status = EXIT_FAILURE;
		printf(" *** error in relxill (%s): %s!\n", func, msg);
	}
}

void get_version_number(char** vstr, int* status){
	if (asprintf(vstr, "%i.%i.%i", version_major, version_minor, version_build) == -1){
		RELXILL_ERROR("failed to get version number",status);
	}
}


/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double* ener, int n_ener, double emin, double emax){
	for (int ii=0; ii<n_ener; ii++){
		ener[ii] = 1.0*ii / (n_ener-1) * ( log(emax) - log(emin)) + log(emin);
		ener[ii] = exp(ener[ii]);
	}
}


