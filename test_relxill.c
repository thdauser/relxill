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

#include "test_relxill.h"

void set_std_param_xillver(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 300.0;  // Ecut
  inp_par[3] = 0.0;    // logxi
  inp_par[4] = 0.;     // redshift
  inp_par[5] = 45.0;   // inclination
  inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_xillverns(double *inp_par) {
  inp_par[0] = 2.0;    // kTB
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 15.0;  // logN
  inp_par[3] = 1.0;    // logxi
  inp_par[4] = 0.;     // redshift
  inp_par[5] = 45.0;   // inclination
  inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_xillverco(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 5.0;    // A_CO
  inp_par[2] = 0.1;   // kTbb
  inp_par[3] = 0.1;   // frac_pl_bb
  inp_par[4] = 300.0;  // Ecut
  inp_par[5] = 0.;     // redshift
  inp_par[6] = 45.0;   // inclination
  inp_par[7] = -1.0;   // refl. frac.
}

void set_std_param_xillver_dens_nthcomp(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 60.0;  // kTe
  inp_par[3] = 15;  // logN
  inp_par[4] = 0.0;    // logxi
  inp_par[5] = 0.;     // redshift
  inp_par[6] = 45.0;   // inclination
  inp_par[7] = -1.0;   // refl. frac.
}

void set_std_param_xillver_nthcomp(double *inp_par) {
  inp_par[0] = 2.1;    // Gamma
  inp_par[1] = 1.0;    // Afe
  inp_par[2] = 60.0;  // kTe
  inp_par[3] = 0.0;    // logxi
  inp_par[4] = 0.;     // redshift
  inp_par[5] = 45.0;   // inclination
  inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_relline(double *inp_par) {
  inp_par[0] = 6.4;
  inp_par[1] = 3.0;
  inp_par[2] = 3.0;
  inp_par[3] = 15.0;
  inp_par[4] = 0.998;
  inp_par[5] = 30.0;
  inp_par[6] = -1.1;
  inp_par[7] = 1000.;
  inp_par[8] = 0.0;
  inp_par[9] = 1.0;
}

void set_std_param_relconv(double *inp_par) {
  inp_par[0] = 3.0;
  inp_par[1] = 3.0;
  inp_par[2] = 15.0;
  inp_par[3] = 0.998;
  inp_par[4] = 30.0;
  inp_par[5] = -1.;
  inp_par[6] = 400.;
  inp_par[7] = 0.0;
}

void set_std_param_relconvlp(double *inp_par) {
  inp_par[0] = 3.0;
  inp_par[1] = 0.998;
  inp_par[2] = 30.0;
  inp_par[3] = -1.;
  inp_par[4] = 400.;
  inp_par[5] = 0.0;
  inp_par[6] = 2.0;
}

void set_std_param_relxill(double *inp_par) {
  inp_par[0] = 3.0;
  inp_par[1] = 3.0;
  inp_par[2] = 15.0;
  inp_par[3] = 0.998;
  inp_par[4] = 60.0;
  inp_par[5] = -1.0;
  inp_par[6] = 400.;
  inp_par[7] = 0.0;   // redshift
  inp_par[8] = 2.1;   // pl Index
  inp_par[9] = 0.0;   // logxi
  inp_par[10] = 1.0;   // Afe
  inp_par[11] = 300.0; // Ecut
  inp_par[12] = 3.0;   // refl_frac
}

void set_std_param_relxill_nthcomp(double *inp_par) {
  inp_par[0] = 3.0;
  inp_par[1] = 3.0;
  inp_par[2] = 15.0;
  inp_par[3] = 0.998;
  inp_par[4] = 60.0;
  inp_par[5] = -1.0;
  inp_par[6] = 400.;
  inp_par[7] = 0.0;   // redshift
  inp_par[8] = 2.1;   // pl Index
  inp_par[9] = 0.0;   // logxi
  inp_par[10] = 1.0;   // Afe
  inp_par[11] = 100.0; // kTe
  inp_par[12] = 0.0;   // refl_frac
}

void set_std_param_relxilldens(double *inp_par) {
  inp_par[0] = 3.0;
  inp_par[1] = 3.0;
  inp_par[2] = 15.0;
  inp_par[3] = 0.998;
  inp_par[4] = 60.0;
  inp_par[5] = -1.0;
  inp_par[6] = 400.;
  inp_par[7] = 0.0;   // redshift
  inp_par[8] = 2.1;   // pl Index
  inp_par[9] = 0.0;   // logxi
  inp_par[10] = 1.0;   // Afe
  inp_par[11] = 15.0; // logN
  inp_par[12] = 3.0;   // refl_frac
}

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

void set_std_param_relxilllp_nthcomp(double *inp_par) {
  inp_par[0] = -1.1;   // height
  inp_par[1] = 0.9798; // a
  inp_par[2] = 60.0;  // incl
  inp_par[3] = -1.0;  // rin
  inp_par[4] = 400.;  // rout
  inp_par[5] = 0.0;    // redshift
  inp_par[6] = 2.1;   // pl Index
  inp_par[7] = 0.0;   // logxi
  inp_par[8] = 1.0;   // Afe
  inp_par[9] = 100.0; // kTe
  inp_par[10] = 3.0;   // refl_frac
  inp_par[11] = 0.0;   // fixReflFrac
}

void set_std_param_relxilllpdens(double *inp_par) {
  inp_par[0] = 3.0;   // height
  inp_par[1] = 0.998; // a
  inp_par[2] = 60.0;  // incl
  inp_par[3] = -1.0;  // rin
  inp_par[4] = 400.;  // rout
  inp_par[5] = 0.0;    // redshift
  inp_par[6] = 2.1;   // pl Index
  inp_par[7] = 0.0;   // logxi
  inp_par[8] = 1.0;   // Afe
  inp_par[9] = 16.0; // logN
  inp_par[10] = -1.0;   // refl_frac
  inp_par[11] = 1.0;   // fixReflFrac
}

void set_std_param_relxilllpdens_nthcomp(double *inp_par) {
  inp_par[0] = 0.998; // a
  inp_par[1] = 60.0;  // incl
  inp_par[2] = -1.0;  // rin
  inp_par[3] = 400.;  // rout
  inp_par[4] = 6.0;   // height
  inp_par[5] = 0.0;   // htop
  inp_par[6] = 0.0;   // beta
  inp_par[7] = 2.1;   // pl Index
  inp_par[8] = 3.1;   // logxi
  inp_par[9] = 1.0;   // Afe
  inp_par[10] = 100.0; // kTe
  inp_par[11] = 15.0; // logN
  inp_par[12] = 3.0;   // refl_frac
  inp_par[13] = 0.0;   // fixReflFrac
  inp_par[14] = 0.0;    // redshift
}

void set_std_param_relline_lp(double *inp_par) {
  inp_par[0] = 1.0;
  inp_par[1] = 3.0;
  inp_par[2] = 0.998;
  inp_par[3] = 30.0;
  inp_par[4] = -1.;
  inp_par[5] = 400.;
  inp_par[6] = 0.0;  // redshift
  inp_par[7] = 0.0;
  inp_par[8] = 2.0;  // gamma
}

void set_std_param_relxilllpion(double *inp_par) {
  inp_par[0] = 6;   // height
  inp_par[1] = 0.9; // a
  inp_par[2] = 30.0;  // incl
  inp_par[3] = -1.0;  // rin
  inp_par[4] = 1000.;  // rout
  inp_par[5] = 0.0;    // redshift
  inp_par[6] = 2.0;   // pl Index
  inp_par[7] = 3.00;   // logxi
  inp_par[8] = 1.0;   // Afe
  inp_par[9] = 300.0; // Ecut
  inp_par[10] = 0.0;     // beta
  inp_par[11] = 0;     // ion_grad_type
  inp_par[12] = 0.0;   // ion_grad_index
  inp_par[13] = 3.0;   // refl_frac
  inp_par[14] = 0.0;   // fixReflFrac
}

void set_std_param_relxilllpion_nthcomp(double *inp_par) {
  inp_par[0] = 6;   // height
  inp_par[1] = 0.9; // a
  inp_par[2] = 30.0;  // incl
  inp_par[3] = -1.0;  // rin
  inp_par[4] = 1000.;  // rout
  inp_par[5] = 0.0;    // redshift
  inp_par[6] = 2.0;   // pl Index
  inp_par[7] = 3.00;   // logxi
  inp_par[8] = 1.0;   // Afe
  inp_par[9] = 100.0; // kTe
  inp_par[10] = 0.0;     // beta
  inp_par[11] = 0;     // ion_grad_type
  inp_par[12] = 0.0;   // ion_grad_index
  inp_par[13] = 3.0;   // refl_frac
  inp_par[14] = 0.0;   // fixReflFrac
}

void set_std_param_relxillns(double *inp_par) {
  inp_par[0] = 3.0;
  inp_par[1] = 3.0;
  inp_par[2] = 15.0;
  inp_par[3] = 0.998;
  inp_par[4] = 60.0;
  inp_par[5] = -1.0;
  inp_par[6] = 400.;
  inp_par[7] = 0.0;   // redshift
  inp_par[8] = 1.0;   // kTbb
  inp_par[9] = 2.0;   // logxi
  inp_par[10] = 1.0;   // Afe
  inp_par[11] = 15.0; // logN
  inp_par[12] = 3.0;   // refl_frac
}

void set_std_param_relxill_bbret(double *inp_par) {
  inp_par[0] = 0.998;
  inp_par[1] = 60.0;
  inp_par[2] = -1.0;
  inp_par[3] = 1000.;
  inp_par[4] = 0.0;   // redshift
  inp_par[5] = 1.0;   // kTbb
  inp_par[6] = 2.0;   // logxi
  inp_par[7] = 1.0;   // Afe
  inp_par[8] = 15.0; // logN
  inp_par[9] = 1.0;   // refl_frac
  inp_par[10] = 1;   // fixReflFrac
  inp_par[11] = 1.4;   // shiftTmaxRRad
}

void set_std_param_relxillco(double *inp_par) {
  inp_par[0] = 3.0;
  inp_par[1] = 3.0;
  inp_par[2] = 15.0;
  inp_par[3] = 0.998;
  inp_par[4] = 60.0;
  inp_par[5] = -1.0;
  inp_par[6] = 1000.;  // Rout
  inp_par[7] = 0.0;   // redshift
  inp_par[8] = 2.0;   // gamma
  inp_par[9] = 10.0;  // A_CO
  inp_par[10] = 5.0;  // kTbb
  inp_par[11] = 0.1;  // frac_pl_bb
  inp_par[12] = 300.0;  // Ecut
  inp_par[13] = 3.0;   // refl_frac
}

xillParam *get_std_param_xillver(int *status) {
  int n_param = NUM_PARAM_XILLVER;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillver(inp_par);
  return init_par_xillver(inp_par, n_param, status);
}

xillParam *get_std_param_xillver_co(int *status) {
  int n_param = NUM_PARAM_XILLVERCO;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillverco(inp_par);
  return init_par_xillver_co(inp_par, n_param, status);
}

xillParam *get_std_param_xillver_ns(int *status) {
  int n_param = NUM_PARAM_XILLVERNS;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillverns(inp_par);
  return init_par_xillver_ns(inp_par, n_param, status);
}

xillParam *get_std_param_xillver_nthcomp(int *status) {
  int n_param = NUM_PARAM_XILLVER_NTHCOMP;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillver_nthcomp(inp_par);
  return init_par_xillver_nthcomp(inp_par,
                                  n_param,
                                  status);  // this is bad design, but kTe and Ecut are identical in the code
}

xillParam *get_std_param_xillver_dens_nthcomp(int *status) {
  int n_param = NUM_PARAM_XILLVERDENS_NTHCOMP;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
  set_std_param_xillver_dens_nthcomp(inp_par);
  return init_par_xillver_dens_nthcomp(inp_par, n_param, status);
}

void get_std_param_relxilllpDCp(relParam **rel_param, xillParam **xill_param, int *status) {
  int n_param = NUM_PARAM_RELXILLLPDENS_NTHCOMP;
  double *inp_par = (double *) malloc(sizeof(double) * n_param);
  CHECK_MALLOC_VOID_STATUS(inp_par, status)

  set_std_param_relxilllpdens_nthcomp(inp_par);
  init_par_relxilllp_dens_nthcomp(rel_param, xill_param, inp_par, n_param, status);

}

relParam *get_std_param_relline(int *status) {
  int n_param = NUM_PARAM_RELLINE;
  double inp_par[NUM_PARAM_RELLINE];
  set_std_param_relline(inp_par);
  return init_par_relline(inp_par, n_param, status);
}

void get_std_param_relxilllp(relParam **p_rel_param, xillParam **p_xill_param, int *status) {
  int n_param = NUM_PARAM_RELXILLLP;
  double inp_par[NUM_PARAM_RELXILLLP];
  set_std_param_relxilllp(inp_par);

  init_par_relxilllp(p_rel_param, p_xill_param, inp_par, n_param, status);
  CHECK_STATUS_VOID(*status);

  assert(is_relxill_model((*p_rel_param)->model_type));
}

/** standard evaluation of the relline model **/
void std_eval_relline(int *status, int n) {

  printf("\n ==> Evaluating RELLINE MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELLINE;
  double inp_par[NUM_PARAM_RELLINE];
  set_std_param_relline(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 300;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.5, 8.0);

  /* call the relline model */
  double photar[n_ener];
  int ii;
  if (n > 2) {
    for (ii = 0; ii < n; ii++) {
      inp_par[4] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      tdrelline(ener, n_ener, photar, inp_par, n_param, status);
    }
  } else {
    inp_par[8] = 0.0;
    tdrelline(ener, n_ener, photar, inp_par, n_param, status);
    inp_par[8] = 1.0;
    tdrelline(ener, n_ener, photar, inp_par, n_param, status);
    inp_par[8] = 2.0;
    tdrelline(ener, n_ener, photar, inp_par, n_param, status);
    inp_par[8] = 0.0;
    tdrelline(ener, n_ener, photar, inp_par, n_param, status);
  }
}

/** standard evaluation of the relline model **/
void std_eval_relconv(int *status, int n) {

  printf("\n ==> Evaluating RELCONV MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELCONV;
  double inp_par[NUM_PARAM_RELCONV];
  set_std_param_relconv(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 2000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.05, 10.0);

  /* call the relline model */
  double photar[n_ener];

  int ii;
  double ener_line = 1.0;
  for (ii = 0; ii < n_ener; ii++) {
    photar[ii] = 0.0;
    if ((ener[ii] < ener_line) && ener[ii + 1] > ener_line) {
      photar[ii] = 1.0; // ener[1]-ener[0];
    }
  }

  // test output
  // save_xillver_spectrum(ener,photar,n_ener,"test_relconv_inp_spectrum.dat");
  tdrelconv(ener, n_ener, photar, inp_par, n_param, status);
  // save_xillver_spectrum(ener,photar,n_ener,"test_relconv_out_spectrum.dat");
}

/** standard evaluation of the relline model **/
void std_eval_relconvlp(int *status, int n) {

  printf("\n ==> Evaluating RELCONV_LP MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELCONVLP;
  double inp_par[NUM_PARAM_RELCONVLP];
  set_std_param_relconvlp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 2000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.05, 10.0);

  /* call the relline model */
  double photar[n_ener];

  int ii;
  double ener_line = 1.0;
  for (ii = 0; ii < n_ener; ii++) {
    photar[ii] = 0.0;
    if ((ener[ii] < ener_line) && ener[ii + 1] > ener_line) {
      photar[ii] = 1.0; // ener[1]-ener[0];
    }
  }

  // test output
  // save_xillver_spectrum(ener,photar,n_ener,"test_relconv_inp_spectrum.dat");
  tdrelconvlp(ener, n_ener, photar, inp_par, n_param, status);
  save_xillver_spectrum(ener, photar, n_ener, "test_relconv_out_spectrum.dat");
}

/** standard evaluation of the relxill model **/
void std_eval_relxill(int *status, int n) {

  printf("\n ==> Evaluating RELXILL MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILL;
  double inp_par[NUM_PARAM_RELXILL];
  set_std_param_relxill(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 3000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];

  // evaluate relxill on a grid point
  inp_par[3] = 0.15;
  tdrelxill(ener, n_ener, photar, inp_par, n_param, status);

  int ii;
  for (ii = 0; ii < n; ii++) {
    if (n > 1) {
      inp_par[3] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      inp_par[7] = 1.0 * ii / (n - 1) * 0.1 - 0.0;
      inp_par[7] = 0.0;
      inp_par[9] = 1.0 * ii / (n - 1) * 0.1;
      tdrelxill(ener, n_ener, photar, inp_par, n_param, status);
    }
    // printf(" testing a=%.3f , lxi=%.2f \n",inp_par[3],inp_par[9]);
    //printf(" testing z=%.3f \n",inp_par[7]);
  }
  tdrelxill(ener, n_ener, photar, inp_par, n_param, status);

}

/** standard evaluation of the relxill model **/
void bugtest_eval_relxill(int *status) {

  printf("\n ==> Evaluating RELXILL MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILL;
  double inp_par[NUM_PARAM_RELXILL];
  set_std_param_relxill(inp_par);
  CHECK_STATUS_VOID(*status);

  int n = 1;

  /* create an energy grid */
  int n_ener = 3000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];

  printf(" TESTING if we can evaluate an Inclination NOT tabulated: \n");

  int ii;
  for (ii = 0; ii < n; ii++) {
    inp_par[4] = 89.0;
    inp_par[5] = -1.0 - 1.0 * ii / (n - 1) * 5;
    tdrelxill(ener, n_ener, photar, inp_par, n_param, status);
  }
  tdrelxill(ener, n_ener, photar, inp_par, n_param, status);

  if (*status == EXIT_FAILURE) {
    printf("  .... no!  (so everything is good)\n");
    *status = EXIT_SUCCESS;
  }
}

/** standard evaluation of the relxill model **/
void std_eval_relxill_nthcomp(int *status, int n) {

  printf("\n ==> Evaluating RELXILL NTHCOMP MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILL;
  double inp_par[NUM_PARAM_RELXILL];
  set_std_param_relxill_nthcomp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 3000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];

  /* int ii;
   for (ii=0; ii<n; ii++){
      if (n>1){
          //			inp_par[3] = 1.0*ii/(n-1)*0.998*2 - 0.998;
          inp_par[7] = 1.0*ii/(n-1)*0.1 - 0.0;
          inp_par[7] = 0.0;
          // inp_par[9] = 1.0*ii/(n-1)*0.1;
      }
      // printf(" testing a=%.3f , lxi=%.2f \n",inp_par[3],inp_par[9]);
      printf(" testing z=%.3f \n",inp_par[7]);
  }*/
  tdrelxill_nthcomp(ener, n_ener, photar, inp_par, n_param, status);
  save_xillver_spectrum(ener, photar, n_ener, "test_relxillCp_spectrum.dat");
}

/** standard evaluation of the relxillD model **/
void std_eval_relxilldens(int *status, int n) {

  printf("\n ==> Evaluating RELXILL_DENS MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILL;
  double inp_par[NUM_PARAM_RELXILL];
  set_std_param_relxilldens(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 4000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];
  int ii;
  for (ii = 0; ii < n; ii++) {
    if (n > 1) {
      inp_par[3] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      inp_par[9] = 1.0 * ii / (n - 1) * 4.7;
    }
    tdrelxilldens(ener, n_ener, photar, inp_par, n_param, status);
  }

}

/** standard evaluation of the relxill model **/
void std_eval_relxilllp(int *status, int n) {

  printf("\n ==> Evaluating RELXILLLP MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLLP;
  double inp_par[NUM_PARAM_RELXILLLP];
  set_std_param_relxilllp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 1000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];
  int ii;

  struct timeval start, end;
  long seconds, useconds;
  gettimeofday(&start, NULL);

  for (ii = 0; ii < n; ii++) {
    if (n > 1) {
      inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      //		inp_par[7] = 1.0*ii/(n-1)*4.7;
      printf(" relxilllp: testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
    }
    tdrelxilllp(ener, n_ener, photar, inp_par, n_param, status);
  }
  double sum = 0.0;
  int jj;
  for (jj = 0; jj < n_ener; jj++) {
    sum += photar[jj];
  }
  printf(" integ flux = %.2e \n", sum);

  gettimeofday(&end, NULL);

  if (n > 1) {
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    double mtime = ((seconds) * 1000. + useconds * 0.001) / ((double) n);

    printf("time per relxilllp evaluation: %.1f milli seconds\n", mtime);
  }

  // test output
  // save_xillver_spectrum(ener,photar,n_ener,"test_relxilllp_spectrum.dat");

}

/** standard evaluation of the relxill model **/
void std_eval_relxilllpion(int *status, int n) {

  printf("\n ==> Evaluating RELXILLLPION MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLLPION;
  double inp_par[NUM_PARAM_RELXILLLPION];
  set_std_param_relxilllpion(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 1000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];

  struct timeval start, end;
  long seconds, useconds;
  gettimeofday(&start, NULL);

  int ii;
  for (ii = 0; ii < n; ii++) {
    if (n > 1) {
      inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      //		inp_par[7] = 1.0*ii/(n-1)*4.7;
      printf(" relxilllpion: testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
    }
    tdrelxilllpion(ener, n_ener, photar, inp_par, n_param, status);
  }

  gettimeofday(&end, NULL);

  if (n > 1) {
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    double mtime = ((seconds) * 1000 + useconds * 0.001) / ((double) n);

    printf("time per relxilllpion evaluation: %.1f milli seconds\n", mtime);
  }

}

/** standard evaluation of the relxill model **/
void std_eval_relxilllp_nthcomp(int *status, int n) {

  printf("\n ==> Evaluating RELXILLLP NTHCOMP MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLLP;
  double inp_par[NUM_PARAM_RELXILLLP];
  set_std_param_relxilllp_nthcomp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 100;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];
  int ii;

  for (ii = 0; ii < n; ii++) {
    if (n > 1) {
      inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      //		inp_par[7] = 1.0*ii/(n-1)*4.7;
      printf(" relxilllp (nthcomp): testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
    }
    tdrelxilllp_nthcomp(ener, n_ener, photar, inp_par, n_param, status);

  }

  // save_xillver_spectrum(ener,photar,n_ener,"test_relxilllp_spectrum.dat");

}

/** standard evaluation of the relxill model **/
void std_eval_relxilllpion_nthcomp(int *status, int n) {

  printf("\n ==> Evaluating RELXILLLP ION NTHCOMP MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLLPION;
  double inp_par[NUM_PARAM_RELXILLLPION];
  set_std_param_relxilllpion_nthcomp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 100;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];
  int ii;

  for (ii = 0; ii < n; ii++) {
    if (n > 1) {
      inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      printf(" relxilllpion (nthcomp): testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
    }
    tdrelxilllpion_nthcomp(ener, n_ener, photar, inp_par, n_param, status);
    printf(" -> %e \n", photar[0]);
  }

}

/** standard evaluation of the relxill model **/
void std_eval_relxilllpdens(int *status, int n) {

  printf("\n ==> Evaluating RELXILLLP_DENS MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLLP;
  double inp_par[NUM_PARAM_RELXILLLP];
  set_std_param_relxilllpdens(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 100;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];
  int ii;

  struct timeval start, end;
  long seconds, useconds;
  gettimeofday(&start, NULL);

  for (ii = 0; ii < n; ii++) {
    if (n > 1) {
      inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
      inp_par[7] = 1.0 * ii / (n - 1) * 4.7;
      printf(" relxilllp: testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
    }
    tdrelxilllpdens(ener, n_ener, photar, inp_par, n_param, status);
  }

  gettimeofday(&end, NULL);

  if (n > 1) {
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    double mtime = ((seconds) * 1000 + useconds * 0.001) / ((double) n);

    printf("time per relxilllp evaluation: %.1f milli seconds\n", mtime);
  }

}

/** standard evaluation of the relline model **/
void std_eval_relline_lp(int *status, int n) {

  printf("\n ==> Evaluating RELLINE LP MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELLINELP;
  double inp_par_lp[n_param];
  set_std_param_relline_lp(inp_par_lp);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 600;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1.5);

  /* call the relline model */
  double photar[n_ener];
  int ii;
  for (ii = 0; ii < n; ii++) {
//		tdrellinelp(ener,n_ener,photar,inp_par_lp,n_param,status);
    inp_par_lp[2] = 0.998;
    tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
    tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
    inp_par_lp[2] = 0.99;
    tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
    inp_par_lp[2] = 0.998;
    tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
    CHECK_STATUS_VOID(*status);
  }

}

void std_eval_xillver(int *status, int n) {

  CHECK_STATUS_VOID(*status);
  printf(" ==> Evaluating Xillver MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_XILLVER;
  double inp_par[n_param];
  set_std_param_xillver(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 2000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.08, 500.0);


  /* call the relline model */
  double photar[n_ener];
  int ii;
  for (ii = 0; ii < n; ii++) {
    tdxillver(ener, n_ener, photar, inp_par, n_param, status);
    CHECK_STATUS_VOID(*status);
  }

}

void std_eval_xillver_nthcomp(int *status, int n) {

  CHECK_STATUS_VOID(*status);
  printf(" ==> Evaluating xillverCp MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_XILLVER_NTHCOMP;
  double inp_par[n_param];
  set_std_param_xillver_nthcomp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 2000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.08, 500.0);

  double photar[n_ener];
  int ii;
  for (ii = 0; ii < n; ii++) {
    tdxillver_nthcomp(ener, n_ener, photar, inp_par, n_param, status);
    CHECK_STATUS_VOID(*status);
  }

}

/** standard evaluation of the relxillNS model **/
void std_eval_relxill_ns(int *status, int n) {

  printf("\n ==> Evaluating RELXILL NS MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLNS;
  double inp_par[NUM_PARAM_RELXILLNS];
  set_std_param_relxillns(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 3000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];

  tdrelxillns(ener, n_ener, photar, inp_par, n_param, status);
}

/** standard evaluation of the relxillNS model **/
void std_eval_relxill_co(int *status, int n) {

  printf("\n ==> Evaluating RELXILL CO MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLCO;
  double inp_par[NUM_PARAM_RELXILLCO];
  set_std_param_relxillco(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 3000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];

  tdrelxillco(ener, n_ener, photar, inp_par, n_param, status);
}

/** standard evaluation of the relxillD model **/
void std_eval_relxilldens_nthcomp(int *status, int n) {

  printf("\n ==> Evaluating relxillDCp MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLDENS_NTHCOMP;
  double inp_par[n_param];
  set_std_param_relxilldens(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 4000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];
  tdrelxilldens_nthcomp(ener, n_ener, photar, inp_par, n_param, status);

}

/** standard evaluation of the relxillD model **/
void std_eval_relxilllpdens_nthcomp(int *status, int n) {

  printf("\n ==> Evaluating relxilllpDCp MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_RELXILLLPDENS_NTHCOMP;
  double inp_par[n_param];
  set_std_param_relxilllpdens_nthcomp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 4000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

  /* call the relline model */
  double photar[n_ener];
  tdrelxilllpdens_nthcomp(ener, n_ener, photar, inp_par, n_param, status);

  printf("  -> %e \n", photar[0]);
  assert(photar[0] > 1e-10);

}

void std_eval_xillver_dens_nthcomp(int *status, int n) {

  CHECK_STATUS_VOID(*status);
  printf(" ==> Evaluating xillverDCp MODEL \n");
  /* set the parameters */
  int n_param = NUM_PARAM_XILLVERDENS_NTHCOMP;
  double inp_par[n_param];
  set_std_param_xillver_dens_nthcomp(inp_par);
  CHECK_STATUS_VOID(*status);

  /* create an energy grid */
  int n_ener = 2000;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.08, 500.0);


  /* call the relline model */
  double photar[n_ener];
  int ii;
  for (ii = 0; ii < n; ii++) {
    tdxillverdens_nthcomp(ener, n_ener, photar, inp_par, n_param, status);
    CHECK_STATUS_VOID(*status);
  }

}

xill_spec *get_std_xill_spec(int *status) {
  xillParam *xill_param = get_std_param_xillver(status);
  xill_spec *xill_spec = get_xillver_spectra(xill_param, status);
  return xill_spec;
}

void get_RelProfileConstEmisZones(rel_spec **p_rel_profile, relParam **p_rel_param, int nzones, int *status) {

  // relline is per default re-normalized, but we need to test the physical normalization here
  putenv("RELLINE_PHYSICAL_NORM=1");

  *p_rel_param = get_std_param_relline(status);
  (*p_rel_param)->num_zones = nzones;
  // set emissivity to constant
  (*p_rel_param)->emis1 = 0.0;
  (*p_rel_param)->emis2 = 0.0;

  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);

  *p_rel_profile = relbase(ener, n_ener, *p_rel_param, NULL, status);
  putenv("RELLINE_PHYSICAL_NORM=0");

}

rel_spec *get_stdRelProfile(int *status) {

  relParam *rel_param = get_std_param_relline(status);

  int n_ener;
  double *ener;
  get_std_relxill_energy_grid(&n_ener, &ener, status);

  return relbase(ener, n_ener, rel_param, NULL, status);
}

void init_std_relXill_spec(rel_spec **rel_profile, double **xill_spec_output, int *status) {

  CHECK_STATUS_VOID(*status);

  *rel_profile = get_stdRelProfile(status);

  double *xill_flux = malloc(sizeof(double) * (*rel_profile)->n_ener);
  CHECK_MALLOC_VOID_STATUS(xill_flux, status);

  xill_spec *xill_spec_table = get_std_xill_spec(status);
  rebin_spectrum((*rel_profile)->ener,
                 xill_flux,
                 (*rel_profile)->n_ener,
                 xill_spec_table->ener,
                 xill_spec_table->flu[0],
                 xill_spec_table->n_ener);
  *xill_spec_output = xill_flux;

}

double calc_FluxInStdBand(const double *flux, double *ener, const int n) {

  const double ELO_REF_RELPROFILE = 0.1;
  const double EHI_REF_RELPROFILE = 1.0;

  return calcSumInEnergyBand(flux, n, ener, ELO_REF_RELPROFILE, EHI_REF_RELPROFILE);
}

double calc_RelatFluxInStdBand(const rel_spec *spec) {

  assert(spec->n_zones == 1);
  return calc_FluxInStdBand(spec->flux[0], spec->ener, spec->n_ener);
}

double calc_XillverFluxInStdBand(const double *xillverSpec, double *ener, int n) {

  return calc_FluxInStdBand(xillverSpec, ener, n);
}

void compareReferenceFlux(double flux, double refFlux, int *status) {

  const double PREC = 1e-6;

  CHECK_STATUS_VOID(*status);

  if (fabs(flux - refFlux) > PREC) {
    RELXILL_ERROR("failed comparing the calculated flux to the reference flux", status);
    printf("  expecting a flux of %e, but calculated %e \n", refFlux, flux);
  }

}

int test_stdEvaluationRelatFlux(const rel_spec *rel_profile) {

  PRINT_RELXILL_TEST_MSG_DEFAULT();

  int status = EXIT_SUCCESS;

  const double ReferenceRelatStdFlux = 8.371512e-01;
  double relatFlux = calc_RelatFluxInStdBand(rel_profile);
  compareReferenceFlux(relatFlux, ReferenceRelatStdFlux, &status);

  print_relxill_test_result(status);
  return status;
}

int test_stdEvaluationXillverFlux(const double *xill_spec, const rel_spec *rel_profile) {

  PRINT_RELXILL_TEST_MSG_DEFAULT();

  int status = EXIT_SUCCESS;

  const double ReferenceXillverStdFlux = 1.802954e+01;
  double xillFlux = calc_XillverFluxInStdBand(xill_spec, rel_profile->ener, rel_profile->n_ener);
  compareReferenceFlux(xillFlux, ReferenceXillverStdFlux, &status);

  print_relxill_test_result(status);
  return status;
}

void test_stdEvaluationFluxes(int *status) {

  CHECK_STATUS_VOID(*status);

  PRINT_RELXILL_TEST_MSG("  : \n");

  rel_spec *rel_profile = NULL;
  double *xill_spec = NULL;
  init_std_relXill_spec(&rel_profile, &xill_spec, status);

  int relStatus = test_stdEvaluationRelatFlux(rel_profile);

  int xillStatus = test_stdEvaluationXillverFlux(xill_spec, rel_profile);

  if ((relStatus != EXIT_SUCCESS) || (xillStatus != EXIT_SUCCESS)) {
    *status = EXIT_FAILURE;
  }
}

