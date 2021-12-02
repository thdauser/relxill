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

    Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef MODELS_H_
#define MODELS_H_

#include "relbase.h"
#include "relutility.h"
#include "xilltable.h"

/**** DEFINES **/
#define MOD_TYPE_RELLINE 1
#define MOD_TYPE_RELLINELP 2

#define MOD_TYPE_RELCONV 11
#define MOD_TYPE_RELCONVLP 12

#define MOD_TYPE_XILLVER 0
#define MOD_TYPE_XILLVER_NTHCOMP 100

#define MOD_TYPE_RELXILL -1
#define MOD_TYPE_RELXILLLP -2

/** density models **/
#define MOD_TYPE_RELXILLDENS -10
#define MOD_TYPE_RELXILLLPDENS -11
#define MOD_TYPE_XILLVERDENS -100

/** ion grad models **/
#define MOD_TYPE_RELXILLLPION -21

/** CO models **/
#define MOD_TYPE_RELXILLCO -200
#define MOD_TYPE_XILLVERCO -210

/** Neutron Star / BB models **/
#define MOD_TYPE_RELXILLNS -30
#define MOD_TYPE_XILLVERNS -101

/** Density Models **/
#define MOD_TYPE_XILLVERDENS_NTHCOMP 1000
#define MOD_TYPE_RELXILLDENS_NTHCOMP -1001
#define MOD_TYPE_RELXILLLPDENS_NTHCOMP -1002


// unpublished models
#define MOD_TYPE_RELXILLBBRET -300
#define MOD_TYPE_RELXILLALPHA -2000
#define MOD_TYPE_RELXILLLPALPHA -2001

/**** FUNCTION DEFINITIONS ****/

void check_parameter_bounds(relParam *param, int *status);
double *shift_energ_spec_1keV(const double *ener, int n_ener, double line_energ, double z, int *status);

/** basic xillver model function **/
void xillver_base(double *ener0, int n_ener0, double *photar, xillParam *param_struct, int *status);

void relline_base(double *ener1keV, double *photar,int n_ener, relParam *param_struct, int *status);


/* get the version number text on the screen (if not already printed before */
void print_version_number(void);

/* get a new relbase parameter structure and initialize it */
relParam *new_relParam(int model_type, int emis_type, int *status);

/* free relbase parameter */
void free_relParam(relParam *);

xillParam *new_xillParam(int model_type, int prim_type, int *status);
void free_xillParam(xillParam *);


#endif /* MODELS_H_ */
