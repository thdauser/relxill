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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef MODELS_H_
#define MODELS_H_

#include "Relbase.h"
#include "ModelParams.h"

extern "C" {
#include "xilltable.h"
#include "relutility.h"
}


/**** FUNCTION DEFINITIONS ****/

void check_parameter_bounds(relParam *param, int *status);

/** basic xillver model function **/
void xillver_base(double *ener0, int n_ener0, double *photar, const ModelParams &inp_params, int *status);

void relline_base(double *ener1keV, double *photar,int n_ener, relParam *param_struct, int *status);

/* get the version number text on the screen (if not already printed before */
void print_version_number(void);

#endif /* MODELS_H_ */
