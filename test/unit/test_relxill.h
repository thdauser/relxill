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
#ifndef RELXILL_TEST_RELXILL_H
#define RELXILL_TEST_RELXILL_H

#include "relbase.h"
#include "relmodels.h"
#include "relutility.h"
#include "xilltable.h"

#include <sys/time.h> // onyl for speed tests

/*
void set_std_param_xillver(double *inp_par);

void set_std_param_xillver_nthcomp(double *inp_par);

void set_std_param_xillverco(double *inp_par);

void set_std_param_xillverns(double *inp_par);

void set_std_param_relline(double *inp_par);

void set_std_param_relxill_bbret(double *inp_par);

xillParam *get_std_param_xillver(int *status);

xillParam *get_std_param_xillver_co(int *status);

xillParam *get_std_param_xillver_ns(int *status);

xillParam *get_std_param_xillver_dens(int *status);

xillParam *get_std_param_xillver_nthcomp(int *status);

xillParam *get_std_param_xillver_dens_nthcomp(int *status);

relParam *get_std_param_relline(int *status);

relParam *get_std_param_rellinelp(int *status);

void get_std_param_relxilllpDCp(relParam **rel_param, xillParam **xill_param, int *status);

 */

xillSpec *get_std_xill_spec(int *status);
rel_spec *get_stdRelProfile(int *status);
void get_RelProfileConstEmisZones(rel_spec **p_rel_profile, relParam **p_rel_param, int nzones, int *status);

void init_std_relXill_spec(rel_spec **rel_profile, double **xill_spec_output, int *status);

void test_stdEvaluationFluxes(int *status);

#endif //RELXILL_TEST_RELXILL_H
