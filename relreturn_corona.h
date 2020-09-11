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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#ifndef RELXILL_RELRETURN_C_RELRETURN_CORONA_H_
#define RELXILL_RELRETURN_C_RELRETURN_CORONA_H_

#include "relbase.h"
#include "common.h"
#include "relmodels.h"

#define NUM_PARAM_RELXILLLPRET 13

void interpolEmisProfile(emisProfile *emisReb, emisProfile *emis0, int *status);

void tdrelxilllpret(const double *ener0,
                    int n_ener0,
                    double *photar,
                    const double *parameter,
                    int n_parameter,
                    int *status);

#endif //RELXILL_RELRETURN_C_RELRETURN_CORONA_H_
