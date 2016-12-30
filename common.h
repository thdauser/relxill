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
#ifndef COMMON_H_
#define COMMON_H_

/*********** DEFINE STATEMENTS *********/

/** define Emissivity Model Type **/
#define EMIS_TYPE_BKN 1
#define EMIS_TYPE_LP 2

/***************************************/

/****** TYPE DEFINITIONS ******/

typedef struct{
	int model_type;
	int emis_type;
	double a;
	double incl;
	double emis1;
	double emis2;
	double rbr;
	double rin;
	double rout;
	double lineE;
	double z;
	double height;
	double gamma;
} relParam;


typedef struct{
	int nr;
	int ng;

	double* re;
	double* gmin;
	double* gmax;
	double* gstar;

	double* emis;

	double*** trff;
	double*** cosne;
} relSysPar;


typedef struct{
	int n_ener;
	int n_zones;
	double* rgrid; // length=n_zones + 1
	double* ener;  // length=n_ener +1
	double** flux;
} rel_spec;

/******************************/


#endif /* COMMON_H_ */
