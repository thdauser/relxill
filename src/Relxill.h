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
#ifndef RELXILL_H_
#define RELXILL_H_



extern "C" {
#include "relmodels.h"
#include "relbase.h"
#include "writeOutfiles.h"
#include "relutility.h"
#include "xillspec.h"
}


void relxill_kernel(double *ener_inp,
                    double *spec_inp,
                    int n_ener_inp,
                    xillParam *xill_param,
                    relParam *rel_param,
                    int *status);

#endif
