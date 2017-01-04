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

/** calculated the emissivity profile, depending on the EMIS_TYPE given **/
void calc_emis_profile(relParam* param, relSysPar* sysPar, int* status){

	assert(sysPar!=NULL);

	if (param->emis_type==EMIS_TYPE_BKN){
		get_emis_bkn(sysPar->emis, sysPar->re, sysPar->nr,
				param->emis1,param->emis2,param->rbr);
	} else if (param->emis_type==EMIS_TYPE_BKN){

	} else {
		RELXILL_ERROR(" calculation of emissivity profile not possible \n",status);
		printf("   -> emis_type=%i not known \n",param->emis_type);
		return;
	}
	return;
}
