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
#include "relutility.h"


/** global parameters, which can be used for several calls of the model */
relParam* cached_rel_param=NULL;
relTable* relline_table=NULL;


/* the relbase function calculating the basic relativistic line shape for a given parameter setup
 * input: ener(n_ener), param
 * output: photar(n_ener)     */
void relbase(const double* ener, const int n_ener, double* photar, const relParam* param, int* status){

	printf(" starting function RELBASE \n");
	// initialize parameter values

	relTable* relline_table=NULL;

	// load tables
	if (relline_table==NULL){
		printf("reading RELLINE table \n");
//		relTable* relline_table=(relTable*) malloc(sizeof(relTable));
		read_relline_table(RELTABLE_FILENAME,&relline_table,status);
	}

	free_relTable(relline_table);
	// get line shape

	// [...]

}


