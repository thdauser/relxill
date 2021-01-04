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
#include "relbase.h"
#include "relmodels.h"


void set_std_param_relxilllp(double *inp_par) {
  inp_par[0] = -1.1;   // height
  inp_par[1] = 0.998; // a
  inp_par[2] = 60.0;  // incl
  inp_par[3] = -1.0;  // rin
  inp_par[4] = 1000.;  // rout
  inp_par[5] = 0.0;    // redshift
  inp_par[6] = 2.1;   // pl Index
  inp_par[7] = 0.0;   // logxi
  inp_par[8] = 1.0;   // Afe
  inp_par[9] = 300.0; // Ecut
  inp_par[10] = 3.0;   // refl_frac
  inp_par[11] = 0.0;   // fixReflFrac
}

void std_eval_relxilllp(int *status) {


  int n_param = NUM_PARAM_RELXILLLP;
  double inp_par[NUM_PARAM_RELXILLLP];
  set_std_param_relxilllp(inp_par);
  CHECK_STATUS_VOID(*status);

  int n_ener = 100;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  double photar[n_ener];
  tdrelxilllp(ener, n_ener, photar, inp_par, n_param, status);

  printf("\n  -> Evaluating RELXILLLP MODEL ... ");
  if (*status == EXIT_SUCCESS) {
    printf("  SUCCESSFUL \n");
  } else {
    printf("  FAILED \n");
  }

}


int main(int argc, char *argv[]) {

  int status = EXIT_SUCCESS;

  if ( (argc >= 2) && (strcmp(argv[1], "version") == 0) ) {
    printf("%s",PROJECT_VER);
  } else {
    std_eval_relxilllp(&status);
  }

  return status;
}
