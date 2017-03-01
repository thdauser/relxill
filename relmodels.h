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

    Copyright 2017 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef MODELS_H_
#define MODELS_H_

#include "relbase.h"
#include "relutility.h"
#include "xilltable.h"

/**** DEFINES **/
#define MOD_TYPE_RELLINE 1
#define NUM_PARAM_RELLINE 9

#define MOD_TYPE_RELLINELP 2
#define NUM_PARAM_RELLINELP 8


#define MOD_TYPE_RELCONV 11
#define NUM_PARAM_RELCONV 8

#define MOD_TYPE_RELCONVLP 12
#define NUM_PARAM_RELCONVLP 7


#define MOD_TYPE_XILLVER 0
#define NUM_PARAM_XILLVER 7

#define MOD_TYPE_RELXILL -1
#define NUM_PARAM_RELXILL 13


/****  TYPE DEFINITIONS ****/


/**** FUNCTION DEFINITIONS ****/
relParam* init_par_relline(const double* inp_par, const int n_parameter, int* status);
relParam* init_par_relline_lp(const double* inp_par, const int n_parameter, int* status);
relParam* init_par_relconv(const double* inp_par, const int n_parameter, int* status);
xillParam* init_par_xillver(const double* inp_par, const int n_parameter, int* status);
void init_par_relxill(relParam** rel_param, xillParam** xill_param, const double* inp_par, const int n_parameter, int* status);



void relline(const double* ener, const int n_ener, double* photar, const double* parameter, const int n_parameter, int* status);

void rellinelp(const double* ener, const int n_ener, double* photar, const double* parameter, const int n_parameter, int* status);

void relxill(const double* ener0, const int n_ener0, double* photar, const double* parameter, const int n_parameter, int* status);

void xillver(const double* ener0, const int n_ener0, double* photar, const double* parameter, const int n_parameter, int* status);

void relconv(const double* ener, const int n_ener, double* photar, const double* parameter, const int n_parameter, int* status);


/* get a new relbase parameter structure and initialize it */
relParam* new_relParam(int model_type, int emis_type, int* status);

/* free relbase parameter */
void free_relParam(relParam*);

xillParam* new_xillParam(int model_type, int prim_type, int* status);
void free_xillParam(xillParam*);

#endif /* MODELS_H_ */
