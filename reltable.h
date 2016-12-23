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
#ifndef RELTABLE_H_
#define RELTABLE_H_

#include "relbase.h"
#include "relutility.h"



/** a single element in the RELLINE table array */
typedef struct{
	float* r;
	float* gmin;
	float* gmax;
	float** trff1;
	float** trff2;
	float** cosne1;
	float** cosne2;
}relDat;

/** the RELLINE table structure */
typedef struct{
	float* a; // spin
	int n_a;

	float* mu0; // inclination
	int n_mu0;

	relDat*** arr; // relline data array

	// dimensions of relline array
	int n_r;
	int n_g;

}relTable;


/* create a new relline table structure */
relTable* new_relTable(int n_a, int n_mu0, int n_r, int n_g, int* status);

/* destroy the relline table structure */
void free_relTable(relTable* tab);

/* routine to read the RELLINE table */
void read_relline_table(char* filename, relTable** tab, int* status);


#endif /* RELTABLE_H_ */
