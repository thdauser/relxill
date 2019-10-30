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

    Copyright 2019 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELXILL_TEST_RELXILL_H
#define RELXILL_TEST_RELXILL_H

#include "relbase.h"
#include "relmodels.h"
#include "relutility.h"
#include "xilltable.h"

#include <sys/time.h> // onyl for speed tests

void set_std_param_xillver(double *inp_par);

void set_std_param_xillverco(double *inp_par);

void set_std_param_relline(double *inp_par);

void set_std_param_relconv(double *inp_par);

void set_std_param_relconvlp(double *inp_par);

void set_std_param_relxill(double *inp_par);

void set_std_param_relxill_nthcomp(double *inp_par);

void set_std_param_relxilldens(double *inp_par);

void set_std_param_relxilllp(double *inp_par);

void set_std_param_relxilllp_nthcomp(double *inp_par);

void set_std_param_relxilllpdens(double *inp_par);

void set_std_param_relline_lp(double *inp_par);

void set_std_param_relxilllpion(double *inp_par);

void set_std_param_relxilllpion_nthcomp(double *inp_par);

void set_std_param_relxillns(double *inp_par);

void set_std_param_relxillco(double *inp_par);

xillParam *get_std_param_xillver(int *status);

xillParam *get_std_param_xillver_co(int *status);

/** standard evaluation of the relline model **/
void std_eval_relline(int *status, int n);

/** standard evaluation of the relline model **/
void std_eval_relconv(int *status, int n);

/** standard evaluation of the relline model **/
void std_eval_relconvlp(int *status, int n);

/** standard evaluation of the relxill model **/
void std_eval_relxill(int *status, int n);

/** standard evaluation of the relxill model **/
void std_eval_relxill_nthcomp(int *status, int n);

/** standard evaluation of the relxillD model **/
void std_eval_relxilldens(int *status, int n);

/** standard evaluation of the relxill model **/
void std_eval_relxilllp(int *status, int n);

/** standard evaluation of the relxill model **/
void std_eval_relxilllpion(int *status, int n);

/** standard evaluation of the relxill model **/
void std_eval_relxilllp_nthcomp(int *status, int n);

/** standard evaluation of the relxill model **/
void std_eval_relxilllpion_nthcomp(int *status, int n);

/** standard evaluation of the relxill model **/
void std_eval_relxilllpdens(int *status, int n);

/** standard evaluation of the relline model **/
void std_eval_relline_lp(int *status, int n);

/** standard evaluation of the relline model **/
void std_eval_xillver(int *status, int n);

/** standard evaluation of the relxillNS model **/
void std_eval_relxill_ns(int *status, int n);

/** standard evaluation of the relxillCO model **/
void std_eval_relxill_co(int *status, int n);

#endif //RELXILL_TEST_RELXILL_H
