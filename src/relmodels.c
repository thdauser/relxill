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

#include "relmodels.h"

int version_number_printed = 0;

static void setNegativeRadiiToRisco(double *r, double a) {
  if (*r < 0) {
    *r = -1.0 * (*r) * kerr_rms(a);
  }
}

static void setNegativeHeightToRplus(double *h, double a) {
  if (*h < 0) {
    *h = -1.0 * (*h) * kerr_rplus(a);
  }
}

void print_version_number(int *status) {
  if (version_number_printed == 0) {
    printf(" *** loading RELXILL model (version %s) *** \n", PROJECT_VER);
    version_number_printed = 1;
  }
}

int warned_rms = 0;
int warned_height = 0;


void check_parameter_bounds(relParam *param, int *status) {

  // first set the Radii to positive value
  setNegativeRadiiToRisco(&(param->rin), param->a);
  setNegativeRadiiToRisco(&(param->rout), param->a);
  setNegativeRadiiToRisco(&(param->rbr), param->a);

  const double rout_max = 1000.0;

  if (param->rout <= param->rin) {
    printf(" *** relxill error : Rin >= Rout not possible, please set the parameters correctly  \n");
    *status = EXIT_FAILURE;
  }

  double rms = kerr_rms(param->a);
  if (param->rin < rms) {
    if (!warned_rms) {
      printf(" *** relxill warning : Rin < ISCO, resetting Rin=ISCO; please set your limits properly \n");
      warned_rms = 1;
    }
    param->rin = rms;
  }

  if (param->a > 0.9982) {
    printf(" *** relxill error : Spin a > 0.9982, model evaluation failed \n");
    *status = EXIT_FAILURE;
    return;
  }

  if (param->a < -1) {
    printf(" *** relxill error : Spin a < -1, model evaluation failed \n");
    *status = EXIT_FAILURE;
    return;
  }

  if (param->incl < 3 * M_PI / 180 || param->incl > 87 * M_PI / 180) {
    printf(" *** relxill error : incl %.3f  is not in the required range between 3-87 deg, model evaluation failed \n",
           param->incl * 180 / M_PI);
    *status = EXIT_FAILURE;
    return;
  }

  if (param->rout <= param->rin) {
    printf(" *** Error : Rout <= Rin, model evaluation failed \n");
    *status = EXIT_FAILURE;
    return;
  }

  if (param->rout > rout_max) {
    printf(
        " *** Error : Rout=%.2e > %.2e Rg, which is the maximal possible value. Make sure to set your limits properly. \n",
        param->rout,
        rout_max);
    printf("             -> resetting Rout=%.2e\n", rout_max);
    param->rout = rout_max;
  }


  /** check rbr values (only applies to BKN emissivity) **/
  if (param->emis_type == EMIS_TYPE_BKN) {
    if (param->rbr < param->rin) {
      printf(" *** warning : Rbr < Rin, resetting Rbr=Rin; please set your limits properly \n");
      param->rbr = param->rin;
    }

    if (param->rbr > param->rout) {
      printf(" *** warning : Rbr > Rout, resetting Rbr=Rout; please set your limits properly \n");
      param->rbr = param->rout;
    }

  }


  /** check velocity values (only applies to LP emissivity) **/
  if (param->emis_type == EMIS_TYPE_LP) {
    if (param->beta < 0) {
      printf(" *** warning (relxill):  beta < 0 is not implemented   (beta=%.3e\n)", param->beta);
      param->beta = 0.0;
    }
    if (param->beta > 0.99) {
      printf(" *** warning (relxill):  velocity has to be within 0 <= beta < 0.99  (beta=%.3e\n)", param->beta);
      param->beta = 0.99;
    }
  }

  /** check height values (only applies to LP emissivity **/
  if (param->emis_type == EMIS_TYPE_LP) {
    setNegativeHeightToRplus(&(param->height), param->a);
    setNegativeHeightToRplus(&(param->htop), param->a);

    double h_fac = 1.1;
    double r_event = kerr_rplus(param->a);
    if ((h_fac * r_event - param->height) > 1e-4) {
      if (!warned_height) {
        printf(" *** Warning : Lamp post source too close to the black hole (h < %.1f r_event) \n", h_fac);
        printf("      Change to negative heights (h <= -%.1f), if you want to fit in units of the Event Horizon \n",
               h_fac);
        printf("      Height= %.3f  ;  r_event=%.3f \n", param->height, r_event);
        printf("      Setting    h =  1.1*r_event  = %.3f \n", r_event * h_fac);
        param->height = r_event * h_fac;
        warned_height = 1;
      }

    }
  }

}

xillParam *init_par_xillver(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVER, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVER);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = inp_par[2];
  param->lxi = inp_par[3];
  param->dens = 15; // logN
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  return param;
}

xillParam *init_par_xillver_nthcomp(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVER_NTHCOMP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVER_NTHCOMP);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = inp_par[2];   // is kTe internally
  param->lxi = inp_par[3];
  param->dens = 15; // logN
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  return param;
}

xillParam *init_par_xillver_ns(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERNS, PRIM_SPEC_BB, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERNS);

  param->kTbb = inp_par[0];
  param->afe = inp_par[1];
  param->ect = 0.0;        // Ecut does not make sense for a BB spectrum
  param->lxi = inp_par[3];
  param->dens = inp_par[2]; // logN
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  // TODO: check parameter bounds here as well
/*	check_parameter_bounds_xillver(param,status);
	CHECK_STATUS_RET(*status,NULL); */

  return param;
}

xillParam *init_par_xillver_co(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERCO, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERCO);

  param->gam = inp_par[0];
  param->afe = inp_par[1]; // this is A_CO here for the xillverCO model
  param->kTbb = inp_par[2]; //
  param->frac_pl_bb = inp_par[3]; //
  param->ect = inp_par[4];
  param->z = inp_par[5];
  param->incl = inp_par[6]; // is given in degrees !!
  param->refl_frac = inp_par[7];

  param->dens = 17.0;
  param->lxi = 0.0;       // interestingly this model does not have an ionization

  // TODO: check parameter bounds here as well
/*	check_parameter_bounds_xillver(param,status);
	CHECK_STATUS_RET(*status,NULL); */

  return param;
}

xillParam *init_par_xillver_dens(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERDENS, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERDENS);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = 300.0;
  param->dens = inp_par[2];
  param->lxi = inp_par[3];
  param->z = inp_par[4];
  param->incl = inp_par[5]; // is given in degrees !!
  param->refl_frac = inp_par[6];

  // TODO: check parameter bounds here as well
/*	check_parameter_bounds_xillver(param,status);
	CHECK_STATUS_RET(*status,NULL); */

  return param;
}

void init_par_relxill(relParam **rel_param,
                      xillParam **xill_param,
                      const double *inp_par,
                      const int n_parameter,
                      int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILL, EMIS_TYPE_BKN, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILL, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILL);

  param->emis1 = inp_par[0];
  param->emis2 = inp_par[1];
  param->rbr = inp_par[2];
  param->a = inp_par[3];
  param->incl = inp_par[4] * M_PI / 180;
  param->rin = inp_par[5];
  param->rout = inp_par[6];
  param->z = inp_par[7];
  xparam->z = inp_par[7];

  xparam->gam = inp_par[8];
  xparam->lxi = inp_par[9];
  xparam->afe = inp_par[10];
  xparam->ect = inp_par[11];
  xparam->dens = 15;

  xparam->refl_frac = inp_par[12];
  xparam->fixReflFrac = 0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

void init_par_relxill_ns(relParam **rel_param,
                         xillParam **xill_param,
                         const double *inp_par,
                         const int n_parameter,
                         int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLNS, EMIS_TYPE_BKN, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLNS, PRIM_SPEC_BB, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLNS);

  param->emis1 = inp_par[0];
  param->emis2 = inp_par[1];
  param->rbr = inp_par[2];
  param->a = inp_par[3];
  param->incl = inp_par[4] * M_PI / 180;
  param->rin = inp_par[5];
  param->rout = inp_par[6];
  param->z = inp_par[7];
  xparam->z = inp_par[7];

  xparam->kTbb = inp_par[8];
  xparam->lxi = inp_par[9];
  xparam->afe = inp_par[10];
  xparam->ect = 0.0; // Ecut does not make sense for BB
  xparam->dens = inp_par[11];

  xparam->refl_frac = inp_par[12];
  xparam->fixReflFrac = 0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

}

void init_par_relxill_co(relParam **rel_param, xillParam **xill_param, const double *inp_par, const int n_parameter,
                         int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLCO, EMIS_TYPE_BKN, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLCO, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLCO);

  param->emis1 = inp_par[0];
  param->emis2 = inp_par[1];
  param->rbr = inp_par[2];
  param->a = inp_par[3];
  param->incl = inp_par[4] * M_PI / 180;
  param->rin = inp_par[5];
  param->rout = inp_par[6];
  param->z = inp_par[7];
  xparam->z = inp_par[7];

  xparam->gam = inp_par[8];
  xparam->afe = inp_par[9]; // this is A_CO here for the xillverCO model
  xparam->kTbb = inp_par[10]; //
  xparam->frac_pl_bb = inp_par[11]; //
  xparam->ect = inp_par[12];
  xparam->refl_frac = inp_par[13];

  xparam->fixReflFrac = 0;

  xparam->dens = 17.;
  xparam->lxi = 0.0;       // interestingly this model does not have an ionization

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

}

void init_par_relxill_nthcomp(relParam **rel_param,
                              xillParam **xill_param,
                              const double *inp_par,
                              const int n_parameter,
                              int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILL, EMIS_TYPE_BKN, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILL, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILL);

  param->emis1 = inp_par[0];
  param->emis2 = inp_par[1];
  param->rbr = inp_par[2];
  param->a = inp_par[3];
  param->incl = inp_par[4] * M_PI / 180;
  param->rin = inp_par[5];
  param->rout = inp_par[6];
  param->z = inp_par[7];
  xparam->z = inp_par[7];

  xparam->gam = inp_par[8];
  xparam->lxi = inp_par[9];
  xparam->afe = inp_par[10];
  xparam->ect = inp_par[11]; // important: this is kte internally
  xparam->dens = 15;

  xparam->refl_frac = inp_par[12];
  xparam->fixReflFrac = 0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

void init_par_relxilldens(relParam **rel_param,
                          xillParam **xill_param,
                          const double *inp_par,
                          const int n_parameter,
                          int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILL, EMIS_TYPE_BKN, status);
  CHECK_STATUS_VOID(*status);

  /** only set the MODEL TYPE RELXILL DENS for xillver **/
  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLDENS, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILL);

  param->emis1 = inp_par[0];
  param->emis2 = inp_par[1];
  param->rbr = inp_par[2];
  param->a = inp_par[3];
  param->incl = inp_par[4] * M_PI / 180;
  param->rin = inp_par[5];
  param->rout = inp_par[6];
  param->z = inp_par[7];
  xparam->z = inp_par[7];

  xparam->gam = inp_par[8];
  xparam->lxi = inp_par[9];
  xparam->afe = inp_par[10];
  xparam->ect = 300.0;
  xparam->dens = inp_par[11];

  xparam->refl_frac = inp_par[12];
  xparam->fixReflFrac = 0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

void init_par_relxilllp(relParam **rel_param,
                        xillParam **xill_param,
                        const double *inp_par,
                        const int n_parameter,
                        int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLLP, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLP, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLP);

  param->height = inp_par[0];
  param->a = inp_par[1];
  param->incl = inp_par[2] * M_PI / 180;
  param->rin = inp_par[3];
  param->rout = inp_par[4];
  param->z = inp_par[5];
  xparam->z = inp_par[5];

  param->gamma = inp_par[6];
  xparam->gam = inp_par[6];
  xparam->lxi = inp_par[7];
  xparam->afe = inp_par[8];
  xparam->ect = inp_par[9];
  xparam->dens = 15.0;

  xparam->refl_frac = inp_par[10];
  xparam->fixReflFrac = (int) (inp_par[11] + 0.5); // make sure there is nor problem with integer conversion

  param->beta = 0.0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

void init_par_relxilllp_nthcomp(relParam **rel_param,
                                xillParam **xill_param,
                                const double *inp_par,
                                const int n_parameter,
                                int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLLP, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLP);

  param->height = inp_par[0];
  param->a = inp_par[1];
  param->incl = inp_par[2] * M_PI / 180;
  param->rin = inp_par[3];
  param->rout = inp_par[4];
  param->z = inp_par[5];
  xparam->z = inp_par[5];

  param->gamma = inp_par[6];
  xparam->gam = inp_par[6];
  xparam->lxi = inp_par[7];
  xparam->afe = inp_par[8];
  xparam->ect = inp_par[9];  // important: this is kte internally
  xparam->dens = 15.0;

  xparam->refl_frac = inp_par[10];
  xparam->fixReflFrac = (int) (inp_par[11] + 0.5); // make sure there is nor problem with integer conversion

  param->beta = 0.0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

void init_par_relxilllp_dens(relParam **rel_param,
                             xillParam **xill_param,
                             const double *inp_par,
                             const int n_parameter,
                             int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLLP, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLPDENS, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLP);

  param->height = inp_par[0];
  param->a = inp_par[1];
  param->incl = inp_par[2] * M_PI / 180;
  param->rin = inp_par[3];
  param->rout = inp_par[4];
  param->z = inp_par[5];
  xparam->z = inp_par[5];

  param->gamma = inp_par[6];
  xparam->gam = inp_par[6];
  xparam->lxi = inp_par[7];
  xparam->afe = inp_par[8];
  xparam->ect = 300.0;
  xparam->dens = inp_par[9];

  xparam->refl_frac = inp_par[10];
  xparam->fixReflFrac = (int) (inp_par[11] + 0.5); // make sure there is nor problem with integer conversion

  param->beta = 0.0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

xillParam *init_par_xillver_dens_nthcomp(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  xillParam *param = new_xillParam(MOD_TYPE_XILLVERDENS_NTHCOMP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_XILLVERDENS_NTHCOMP);

  param->gam = inp_par[0];
  param->afe = inp_par[1];
  param->ect = inp_par[2];
  param->dens = inp_par[3];
  param->lxi = inp_par[4];
  param->z = inp_par[5];
  param->incl = inp_par[6]; // is given in degrees !!
  param->refl_frac = inp_par[7];

  // TODO: check parameter bounds here as well
/*	check_parameter_bounds_xillver(param,status);
	CHECK_STATUS_RET(*status,NULL); */

  return param;
}

void init_par_relxilldens_nthcomp(relParam **rel_param,
                                  xillParam **xill_param,
                                  const double *inp_par,
                                  const int n_parameter,
                                  int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILL, EMIS_TYPE_BKN, status);
  CHECK_STATUS_VOID(*status);

  /** only set the MODEL TYPE RELXILL DENS for xillver **/
  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLDENS_NTHCOMP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLDENS_NTHCOMP);

  param->emis1 = inp_par[0];
  param->emis2 = inp_par[1];
  param->rbr = inp_par[2];
  param->a = inp_par[3];
  param->incl = inp_par[4] * M_PI / 180;
  param->rin = inp_par[5];
  param->rout = inp_par[6];
  param->z = inp_par[7];
  xparam->z = inp_par[7];

  xparam->gam = inp_par[8];
  xparam->lxi = inp_par[9];
  xparam->afe = inp_par[10];
  xparam->ect = inp_par[11];
  xparam->dens = inp_par[12];

  xparam->refl_frac = inp_par[13];
  xparam->fixReflFrac = 0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

}

void init_par_relxilllp_dens_nthcomp(relParam **rel_param,
                                     xillParam **xill_param,
                                     const double *inp_par,
                                     const int n_parameter,
                                     int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLLP, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLPDENS_NTHCOMP, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLPDENS_NTHCOMP);

  param->a = inp_par[0];
  param->incl = inp_par[1] * M_PI / 180;
  param->rin = inp_par[2];
  param->rout = inp_par[3];
  param->height = inp_par[4];
  param->htop = inp_par[5];
  param->beta = inp_par[6];

  param->gamma = inp_par[7];
  xparam->gam = param->gamma;
  xparam->lxi = inp_par[8];
  xparam->afe = inp_par[9];
  xparam->ect = inp_par[10];
  xparam->dens = inp_par[11];

  xparam->refl_frac = inp_par[12];
  xparam->fixReflFrac = (int) (inp_par[13] + 0.5); // make sure there is no problem with integer conversion

  param->z = inp_par[14];
  xparam->z = param->z;

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

}

relParam *init_par_relline(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELLINE, EMIS_TYPE_BKN, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_RELLINE);

  param->lineE = inp_par[0];
  param->emis1 = inp_par[1];
  param->emis2 = inp_par[2];
  param->rbr = inp_par[3];
  param->a = inp_par[4];
  param->incl = inp_par[5] * M_PI / 180;
  param->rin = inp_par[6];
  param->rout = inp_par[7];
  param->z = inp_par[8];
  param->limb = (int) (inp_par[9] + 0.5);

  check_parameter_bounds(param, status);
  CHECK_STATUS_RET(*status, NULL);

  return param;
}

relParam *init_par_relconv(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELCONV, EMIS_TYPE_BKN, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_RELCONV);

  param->lineE = 1.0;
  param->emis1 = inp_par[0];
  param->emis2 = inp_par[1];
  param->rbr = inp_par[2];
  param->a = inp_par[3];
  param->incl = inp_par[4] * M_PI / 180;
  param->rin = inp_par[5];
  param->rout = inp_par[6];
  param->z = 0.0;
  param->limb = (int) (inp_par[7] + 0.5);

  check_parameter_bounds(param, status);
  CHECK_STATUS_RET(*status, NULL);

  return param;
}

relParam *init_par_relline_lp(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELLINELP, EMIS_TYPE_LP, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_RELLINELP);

  param->lineE = inp_par[0];
  param->height = inp_par[1];
  param->a = inp_par[2];
  param->incl = inp_par[3] * M_PI / 180;
  param->rin = inp_par[4];
  param->rout = inp_par[5];
  param->z = inp_par[6];
  param->limb = (int) (inp_par[7] + 0.5);
  param->gamma = inp_par[8];

  param->beta = 0.0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_RET(*status, NULL);

  return param;
}

relParam *init_par_relconv_lp(const double *inp_par, const int n_parameter, int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELCONVLP, EMIS_TYPE_LP, status);
  CHECK_STATUS_RET(*status, NULL);

  assert(n_parameter == NUM_PARAM_RELCONVLP);

  param->lineE = 1.0;
  param->height = inp_par[0];
  param->a = inp_par[1];
  param->incl = inp_par[2] * M_PI / 180;
  param->rin = inp_par[3];
  param->rout = inp_par[4];
  param->z = 0.0;
  param->limb = (int) (inp_par[5] + 0.5);
  param->gamma = inp_par[6];

  param->beta = 0.0;

  check_parameter_bounds(param, status);
  CHECK_STATUS_RET(*status, NULL);

  return param;
}

void init_par_relxilllpion(relParam **rel_param,
                           xillParam **xill_param,
                           const double *inp_par,
                           const int n_parameter,
                           int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLLPION, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLPION, PRIM_SPEC_ECUT, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLPION);

  param->height = inp_par[0];
  param->a = inp_par[1];
  param->incl = inp_par[2] * M_PI / 180;
  param->rin = inp_par[3];
  param->rout = inp_par[4];
  param->z = inp_par[5];
  xparam->z = inp_par[5];

  param->gamma = inp_par[6];
  xparam->gam = inp_par[6];
  xparam->lxi = inp_par[7];
  xparam->afe = inp_par[8];
  xparam->ect = inp_par[9];

  xparam->dens = 15.0;
  xparam->ion_grad_type = (int) (inp_par[11] + 0.5);  // make sure there is no problem with integer conversion
  xparam->ion_grad_index = inp_par[12];

  xparam->refl_frac = inp_par[13];
  xparam->fixReflFrac = (int) (inp_par[14] + 0.5); // make sure there is no problem with integer conversion

  param->beta = inp_par[10];

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

void init_par_relxilllpion_nthcomp(relParam **rel_param,
                                   xillParam **xill_param,
                                   const double *inp_par,
                                   const int n_parameter,
                                   int *status) {

  // fill in parameters
  relParam *param = new_relParam(MOD_TYPE_RELXILLLPION, EMIS_TYPE_LP, status);
  CHECK_STATUS_VOID(*status);

  xillParam *xparam = new_xillParam(MOD_TYPE_RELXILLLPION, PRIM_SPEC_NTHCOMP, status);
  CHECK_STATUS_VOID(*status);

  assert(n_parameter == NUM_PARAM_RELXILLLPION);

  param->height = inp_par[0];
  param->a = inp_par[1];
  param->incl = inp_par[2] * M_PI / 180;
  param->rin = inp_par[3];
  param->rout = inp_par[4];
  param->z = inp_par[5];
  xparam->z = inp_par[5];

  param->gamma = inp_par[6];
  xparam->gam = inp_par[6];
  xparam->lxi = inp_par[7];
  xparam->afe = inp_par[8];
  xparam->ect = inp_par[9];  // treated as kTe internally (!)

  xparam->dens = 15.0;
  xparam->ion_grad_type = (int) (inp_par[11] + 0.5);  // make sure there is no problem with integer conversion
  xparam->ion_grad_index = inp_par[12];

  xparam->refl_frac = inp_par[13];
  xparam->fixReflFrac = (int) (inp_par[14] + 0.5); // make sure there is no problem with integer conversion

  param->beta = inp_par[10];

  check_parameter_bounds(param, status);
  CHECK_STATUS_VOID(*status);

  *rel_param = param;
  *xill_param = xparam;

  return;
}

/** shift the spectrum such that we can calculate the line for 1 keV **/
double *shift_energ_spec_1keV(const double *ener, const int n_ener, double line_energ, double z, int *status) {

  double *ener1keV = (double *) malloc((n_ener + 1) * sizeof(double));
  CHECK_MALLOC_RET_STATUS(ener1keV, status, NULL)

  int ii;
  for (ii = 0; ii <= n_ener; ii++) {
    ener1keV[ii] = ener[ii] * (1 + z) / line_energ;
  }
  return ener1keV;
}

/** RELXILL MODEL FUNCTION **/
void tdrelxill(const double *ener0,
               const int n_ener0,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxill(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  // int n_ener = (int) n_ener0;
  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_kernel(ener, photar, n_ener0, xill_param, rel_param, status);

  free(ener);
  free_xillParam(xill_param);
  free_relParam(rel_param);

}

/** RELXILL MODEL FUNCTION for the xillverNS table (with BB instead of Power Law) **/
void tdrelxillns(const double *ener0,
                 const int n_ener0,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxill_ns(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_kernel(ener, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);

  free(ener);
  free_xillParam(xill_param);
  free_relParam(rel_param);

}

void tdrelxillco(const double *ener0, const int n_ener0, double *photar, const double *parameter,
                 const int n_parameter, int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxill_co(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  // int n_ener = (int) n_ener0;
  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_kernel(ener, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);

  free(ener);
  free_xillParam(xill_param);
  free_relParam(rel_param);

}

/** RELXILL DENS MODEL FUNCTION **/
void tdrelxilldens(const double *ener0,
                   const int n_ener0,
                   double *photar,
                   const double *parameter,
                   const int n_parameter,
                   int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilldens(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_kernel(ener, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);

  free(ener);
  free_xillParam(xill_param);
  free_relParam(rel_param);

}

/** RELXILL NTHCOMP MODEL FUNCTION **/
void tdrelxill_nthcomp(const double *ener0,
                       const int n_ener0,
                       double *photar,
                       const double *parameter,
                       const int n_parameter,
                       int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxill_nthcomp(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  // int n_ener = (int) n_ener0;
  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_kernel(ener, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);

  free(ener);
  free_xillParam(xill_param);
  free_relParam(rel_param);

}

/** XSPEC RELXILLLP MODEL FUNCTION **/
void tdrelxilllp(const double *ener0,
                 const int n_ener0,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status) {

  CHECK_STATUS_VOID(*status);

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilllp(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  relxill_kernel(ener_shifted, photar, n_ener0, xill_param, rel_param, status);

  free_xillParam(xill_param);
  free_relParam(rel_param);
  free(ener_shifted);

}

/** XSPEC RELXILLLP NTHCOMP MODEL FUNCTION **/
void tdrelxilllp_nthcomp(const double *ener0,
                         const int n_ener0,
                         double *photar,
                         const double *parameter,
                         const int n_parameter,
                         int *status) {

  CHECK_STATUS_VOID(*status);

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilllp_nthcomp(&rel_param, &xill_param, parameter, n_parameter, status);

  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  relxill_kernel(ener_shifted, photar, n_ener0, xill_param, rel_param, status);

  free_xillParam(xill_param);
  free_relParam(rel_param);
  free(ener_shifted);

}

/** XSPEC RELXILLLP MODEL FUNCTION **/
void tdrelxilllpdens(const double *ener0,
                     const int n_ener0,
                     double *photar,
                     const double *parameter,
                     const int n_parameter,
                     int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilllp_dens(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  relxill_kernel(ener_shifted, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);


  free_xillParam(xill_param);
  free_relParam(rel_param);
  free(ener_shifted);

}

/** RELXILL DENS MODEL FUNCTION **/
void tdrelxilldens_nthcomp(const double *ener0,
                           const int n_ener0,
                           double *photar,
                           const double *parameter,
                           const int n_parameter,
                           int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilldens_nthcomp(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);

  relxill_kernel(ener, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);

  free(ener);
  free_xillParam(xill_param);
  free_relParam(rel_param);

}

/** XSPEC RELXILLLP MODEL FUNCTION **/
void tdrelxilllpdens_nthcomp(const double *ener0,
                             const int n_ener0,
                             double *photar,
                             const double *parameter,
                             const int n_parameter,
                             int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilllp_dens_nthcomp(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);


  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  relxill_kernel(ener_shifted, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);


  free_xillParam(xill_param);
  free_relParam(rel_param);
  free(ener_shifted);

}

/** XSPEC XILLVER DENS MODEL FUNCTION **/
void tdxillverdens_nthcomp(const double *ener0,
                           const int n_ener0,
                           double *photar,
                           const double *parameter,
                           const int n_parameter,
                           int *status) {

  xillParam *param_struct = init_par_xillver_dens_nthcomp(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  xillver_base(ener0, n_ener0, photar, param_struct, status);
}

/** XSPEC XILLVER MODEL FUNCTION **/
void tdxillver(const double *ener0,
               const int n_ener0,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status) {

  xillParam *param_struct = init_par_xillver(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  xillver_base(ener0, n_ener0, photar, param_struct, status);
}

/** XSPEC XILLVER NTHCOMP MODEL FUNCTION **/
void tdxillver_nthcomp(const double *ener0,
                       const int n_ener0,
                       double *photar,
                       const double *parameter,
                       const int n_parameter,
                       int *status) {

  xillParam *param_struct = init_par_xillver_nthcomp(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  xillver_base(ener0, n_ener0, photar, param_struct, status);
}

/** XSPEC XILLVER DENS MODEL FUNCTION **/
void tdxillverdens(const double *ener0,
                   const int n_ener0,
                   double *photar,
                   const double *parameter,
                   const int n_parameter,
                   int *status) {

  xillParam *param_struct = init_par_xillver_dens(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  xillver_base(ener0, n_ener0, photar, param_struct, status);
}

/** XSPEC XILLVER NS MODEL FUNCTION **/
void tdxillverns(const double *ener0,
                 const int n_ener0,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status) {

  xillParam *param_struct = init_par_xillver_ns(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  xillver_base(ener0, n_ener0, photar, param_struct, status);
}

/** XSPEC XILLVER CO MODEL FUNCTION **/
void tdxillverco(const double *ener0, const int n_ener0, double *photar,
                 const double *parameter, const int n_parameter, int *status) {

  xillParam *param_struct = init_par_xillver_co(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  xillver_base(ener0, n_ener0, photar, param_struct, status);
}

/** BASIC XILLVER MODEL FUNCTION **/
void xillver_base(const double *ener0, const int n_ener0, double *photar, xillParam *param_struct, int *status) {


  // call the function which calculates the xillver spectrum
  xillSpec *spec = get_xillver_spectra(param_struct, status);
  CHECK_STATUS_VOID(*status);

  // =4= rebin to the input grid
  assert(spec->n_incl == 1); // make sure there is only one spectrum given (for the chosen inclination)

  /** add the dependence on incl, assuming a semi-infinite slab **/
  norm_xillver_spec(spec, param_struct->incl);

  double *ener = (double *) ener0;
  int n_ener = (int) n_ener0;
  double flux[n_ener];

  rebin_spectrum(ener, flux, n_ener, spec->ener, spec->flu[0], spec->n_ener);

  add_primary_component(ener, n_ener, flux, NULL, param_struct, status);

  double *ener_shifted = shift_energ_spec_1keV(ener, n_ener, 1.0, param_struct->z, status);

  rebin_spectrum(ener_shifted, photar, n_ener, ener, flux, n_ener);

  /** free **/
  free_xillParam(param_struct);
  free_xill_spec(spec);
  free(ener_shifted);
}

void relline_base(double *ener1keV, double *photar, const int n_ener, relParam *param_struct, int *status) {

  rel_spec *spec = relbase(ener1keV, n_ener, param_struct, NULL, status);

  for (int ii = 0; ii < n_ener; ii++) {
    photar[ii] = spec->flux[0][ii];
  }
}

/** XSPEC RELLINE MODEL FUNCTION **/
void tdrelline(const double *ener,
               const int n_ener,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status) {

  relParam *param_struct = init_par_relline(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  // shift the spectrum such that we can calculate the line for 1 keV
  double *ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z, status);
  CHECK_STATUS_VOID(*status);

  // call the function which calculates the line (assumes a line at 1keV!)
  relline_base(ener1keV, photar, n_ener, param_struct, status);

  free_relParam(param_struct);
  free(ener1keV);
  CHECK_STATUS_VOID(*status);

}

/** XSPEC RELLINELP MODEL FUNCTION **/
void tdrellinelp(const double *ener,
                 const int n_ener,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status) {

  relParam *param_struct = init_par_relline_lp(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  // shift the spectrum such that we can calculate the line for 1 keV
  double *ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z, status);
  CHECK_STATUS_VOID(*status);

  // call the function which calculates the line (assumes a line at 1keV!)
  rel_spec *spec = relbase(ener1keV, n_ener, param_struct, NULL, status);

  free_relParam(param_struct);
  free(ener1keV);

  CHECK_STATUS_VOID(*status);

  // now we certainly have a valid return value
  assert(spec->n_zones == 1);
  int ii;
  for (ii = 0; ii < n_ener; ii++) {
    photar[ii] = spec->flux[0][ii];
  }

}

/** XSPEC RELLINE MODEL FUNCTION **/
void tdrelconv(const double *ener,
               const int n_ener,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status) {

  relParam *param_struct = init_par_relconv(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  // shift the spectrum such that we can calculate the line for 1 keV
  double *ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z, status);
  CHECK_STATUS_VOID(*status);

  // call the function which calculates the line (assumes a line at 1keV!)
  relconv_kernel(ener1keV, photar, n_ener, param_struct, status);
  CHECK_STATUS_VOID(*status);

  free_relParam(param_struct);
  free(ener1keV);
}

/** XSPEC RELCONV LP MODEL FUNCTION **/
void tdrelconvlp(const double *ener,
                 const int n_ener,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status) {

  relParam *param_struct = init_par_relconv_lp(parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  // shift the spectrum such that we can calculate the line for 1 keV
  double *ener1keV = shift_energ_spec_1keV(ener, n_ener, param_struct->lineE, param_struct->z, status);
  CHECK_STATUS_VOID(*status);

  // call the function which calculates the line (assumes a line at 1keV!)
  relconv_kernel(ener1keV, photar, n_ener, param_struct, status);
  CHECK_STATUS_VOID(*status);

  free_relParam(param_struct);
  free(ener1keV);
}

/** XSPEC RELXILLLPION MODEL FUNCTION **/
void tdrelxilllpion(const double *ener0,
                    const int n_ener0,
                    double *photar,
                    const double *parameter,
                    const int n_parameter,
                    int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilllpion(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);


  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  relxill_kernel(ener_shifted, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);


  free_xillParam(xill_param);
  free_relParam(rel_param);
  free(ener_shifted);

}

/** XSPEC RELXILLLPIONCP MODEL FUNCTION **/
void tdrelxilllpion_nthcomp(const double *ener0,
                            const int n_ener0,
                            double *photar,
                            const double *parameter,
                            const int n_parameter,
                            int *status) {

  xillParam *xill_param = NULL;
  relParam *rel_param = NULL;

  init_par_relxilllpion_nthcomp(&rel_param, &xill_param, parameter, n_parameter, status);
  CHECK_STATUS_VOID(*status);

  double *ener_shifted = shift_energ_spec_1keV(ener0, n_ener0, 1.0, rel_param->z, status);
  relxill_kernel(ener_shifted, photar, n_ener0, xill_param, rel_param, status);
  CHECK_STATUS_VOID(*status);


  free_xillParam(xill_param);
  free_relParam(rel_param);
  free(ener_shifted);

}

/* get a new relbase parameter structure and initialize it */
relParam *new_relParam(int model_type, int emis_type, int *status) {
  relParam *param = (relParam *) malloc(sizeof(relParam));
  if (param == NULL) {
    RELXILL_ERROR("memory allocation failed", status);
    return NULL;
  }
  param->model_type = model_type;
  param->emis_type = emis_type;

  param->a = PARAM_DEFAULT;
  param->incl = PARAM_DEFAULT;
  param->emis1 = PARAM_DEFAULT;
  param->emis2 = PARAM_DEFAULT;
  param->rbr = PARAM_DEFAULT;
  param->rin = PARAM_DEFAULT;
  param->rout = PARAM_DEFAULT;
  param->lineE = PARAM_DEFAULT;
  param->z = PARAM_DEFAULT;
  param->height = 0.0;
  param->gamma = PARAM_DEFAULT;
  param->beta = 0.0; // special case, in order to prevent strange results
  param->height = PARAM_DEFAULT;
  param->htop = 0.0;
  param->limb = 0;
  param->return_rad = 0;

  // this is set by the environment variable "RELLINE_PHYSICAL_NORM"
  param->do_renorm_relline = do_renorm_model(param);

  // set depending on model/emis type and ENV "RELXILL_NUM_RZONES"
  //  -> note as this is onl for relat. models, in case of an ion gradient this needs to be updated
  param->num_zones = get_num_zones(param->model_type, param->emis_type, ION_GRAD_TYPE_CONST);

  return param;
}

/* free relbase parameter */
void free_relParam(relParam *param) {
  free(param);
}

/* get a new relbase parameter structure and initialize it */
xillParam *new_xillParam(int model_type, int prim_type, int *status) {
  xillParam *param = (xillParam *) malloc(sizeof(xillParam));
  if (param == NULL) {
    RELXILL_ERROR("memory allocation failed", status);
    return NULL;
  }
  param->model_type = model_type;
  param->prim_type = prim_type;

  param->gam = PARAM_DEFAULT;
  param->afe = PARAM_DEFAULT;
  param->lxi = PARAM_DEFAULT;
  param->ect = PARAM_DEFAULT;
  param->incl = PARAM_DEFAULT;
  param->z = PARAM_DEFAULT;
  param->refl_frac = PARAM_DEFAULT;
  param->fixReflFrac = -1;
  param->dens = 15;  // the standard value for every table, given in "log-units"
  param->frac_pl_bb = PARAM_DEFAULT;
  param->kTbb = PARAM_DEFAULT;
  param->ion_grad_type = ION_GRAD_TYPE_CONST; // no ion grad
  param->ion_grad_index = PARAM_DEFAULT;

  return param;
}

/* free relbase parameter */
void free_xillParam(xillParam *param) {
  free(param);
}




/***************************************************/
/**** XSPEC MODEL WRAPPERS *************************/
/***************************************************/

//
//// XSPEC RELXILL MODEL FUNCTION
//void lmodrelxill(const double *ener0,
//                 const int n_ener0,
//                 const double *parameter,
//                 int ifl,
//                 double *photar,
//                 double *photer,
//                 const char *init) {
//
//  const int n_parameter = 13;
//  int status = EXIT_SUCCESS;
//
//  tdrelxill(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxill model failed", &status);
//}
//
//// XSPEC RELXILL NS MODEL FUNCTION
//void lmodrelxillns(const double *ener0,
//                   const int n_ener0,
//                   const double *parameter,
//                   int ifl,
//                   double *photar,
//                   double *photer,
//                   const char *init) {
//
//  const int n_parameter = 13;
//  int status = EXIT_SUCCESS;
//  tdrelxillns(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxill_ns model failed", &status);
//}
//
///** XSPEC RELXILL DENS MODEL FUNCTION **/
//void lmodrelxilldens(const double *ener0,
//                     const int n_ener0,
//                     const double *parameter,
//                     int ifl,
//                     double *photar,
//                     double *photer,
//                     const char *init) {
//
//  const int n_parameter = 13;
//  int status = EXIT_SUCCESS;
//  tdrelxilldens(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxill_dens model failed", &status);
//}
//
///** XSPEC RELXILL MODEL FUNCTION **/
//void lmodrelxillnthcomp(const double *ener0,
//                        const int n_ener0,
//                        const double *parameter,
//                        int ifl,
//                        double *photar,
//                        double *photer,
//                        const char *init) {
//
//  const int n_parameter = 13;
//  int status = EXIT_SUCCESS;
//  tdrelxill_nthcomp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxill_nthcomp model failed", &status);
//}
//
///** XSPEC RELXILLLP MODEL FUNCTION **/
//void lmodrelxilllp(const double *ener0,
//                   const int n_ener0,
//                   const double *parameter,
//                   int ifl,
//                   double *photar,
//                   double *photer,
//                   const char *init) {
//
//  const int n_parameter = 12;
//  int status = EXIT_SUCCESS;
//  tdrelxilllp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating rellinelp model failed", &status);
//}
//
///** XSPEC RELXILLLPDENS MODEL FUNCTION **/
//void lmodrelxilllpdens(const double *ener0,
//                       const int n_ener0,
//                       const double *parameter,
//                       int ifl,
//                       double *photar,
//                       double *photer,
//                       const char *init) {
//
//  const int n_parameter = 12;
//  int status = EXIT_SUCCESS;
//  tdrelxilllpdens(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating rellinelp model failed", &status);
//}
//
///** XSPEC RELXILLLP MODEL FUNCTION **/
//void lmodrelxilllpnthcomp(const double *ener0,
//                          const int n_ener0,
//                          const double *parameter,
//                          int ifl,
//                          double *photar,
//                          double *photer,
//                          const char *init) {
//
//  const int n_parameter = 12;
//  int status = EXIT_SUCCESS;
//  tdrelxilllp_nthcomp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating rellinelp model failed", &status);
//}
//
///** XSPEC XILLVER MODEL FUNCTION **/
//void lmodxillver(const double *ener0,
//                 const int n_ener0,
//                 const double *parameter,
//                 int ifl,
//                 double *photar,
//                 double *photer,
//                 const char *init) {
//
//  const int n_parameter = 7;
//  int status = EXIT_SUCCESS;
//  tdxillver(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating xillver model failed", &status);
//}
//
///** XSPEC XILLVER MODEL FUNCTION **/
//void lmodxillverdens(const double *ener0,
//                     const int n_ener0,
//                     const double *parameter,
//                     int ifl,
//                     double *photar,
//                     double *photer,
//                     const char *init) {
//
//  const int n_parameter = 7;
//  int status = EXIT_SUCCESS;
//  tdxillverdens(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating xillver model failed", &status);
//}
//
///** XSPEC XILLVER NS MODEL FUNCTION **/
//void lmodxillverns(const double *ener0,
//                   const int n_ener0,
//                   const double *parameter,
//                   int ifl,
//                   double *photar,
//                   double *photer,
//                   const char *init) {
//
//  const int n_parameter = 7;
//  int status = EXIT_SUCCESS;
//  tdxillverns(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating xillverNS model failed", &status);
//}
//
//void lmodxillverco(const double *ener0, const int n_ener0, const double *parameter, int ifl,
//                   double *photar, double *photer, const char *init) {
//
//  const int n_parameter = 8;
//  int status = EXIT_SUCCESS;
//  tdxillverco(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating xillverCO model failed", &status);
//}
//
///** XSPEC RELXILL CO MODEL FUNCTION **/
//void
//lmodrelxillco(const double *ener0, const int n_ener0, const double *parameter, int ifl, double *photar, double *photer,
//              const char *init) {
//
//  const int n_parameter = 14;
//  int status = EXIT_SUCCESS;
//  tdrelxillco(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxillCO model failed", &status);
//}
//
///** XSPEC XILLVER MODEL FUNCTION **/
//void lmodxillvernthcomp(const double *ener0,
//                        const int n_ener0,
//                        const double *parameter,
//                        int ifl,
//                        double *photar,
//                        double *photer,
//                        const char *init) {
//
//  const int n_parameter = 7;
//  int status = EXIT_SUCCESS;
//  tdxillver_nthcomp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating xillver model failed", &status);
//}
//
///** XSPEC RELLINE MODEL FUNCTION **/
//void lmodrelline(const double *ener0,
//                 const int n_ener0,
//                 const double *parameter,
//                 int ifl,
//                 double *photar,
//                 double *photer,
//                 const char *init) {
//
//  const int n_parameter = 10;
//  int status = EXIT_SUCCESS;
//  tdrelline(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relline model failed", &status);
//}
//
///** XSPEC RELLINELP MODEL FUNCTION **/
//void lmodrellinelp(const double *ener0,
//                   const int n_ener0,
//                   const double *parameter,
//                   int ifl,
//                   double *photar,
//                   double *photer,
//                   const char *init) {
//
//  const int n_parameter = 9;
//  int status = EXIT_SUCCESS;
//  tdrellinelp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating rellinelp model failed", &status);
//}
//
///** XSPEC RELCONV MODEL FUNCTION **/
//void lmodrelconv(const double *ener0,
//                 const int n_ener0,
//                 const double *parameter,
//                 int ifl,
//                 double *photar,
//                 double *photer,
//                 const char *init) {
//
//  const int n_parameter = 8;
//  int status = EXIT_SUCCESS;
//  tdrelconv(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relconv model failed", &status);
//}
//
///** XSPEC RELCONV LP MODEL FUNCTION **/
//void lmodrelconvlp(const double *ener0,
//                   const int n_ener0,
//                   const double *parameter,
//                   int ifl,
//                   double *photar,
//                   double *photer,
//                   const char *init) {
//
//  const int n_parameter = 7;
//  int status = EXIT_SUCCESS;
//  tdrelconvlp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relconv model failed", &status);
//}
//
///** XSPEC RELXILLLP MODEL FUNCTION **/
//void lmodrelxilllpion(const double *ener0,
//                      const int n_ener0,
//                      const double *parameter,
//                      int ifl,
//                      double *photar,
//                      double *photer,
//                      const char *init) {
//
//  const int n_parameter = 15;
//  int status = EXIT_SUCCESS;
//  tdrelxilllpion(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxilllpion model failed", &status);
//}
//
///** XSPEC RELXILLLP MODEL FUNCTION **/
//void lmodrelxilllpionnthcomp(const double *ener0,
//                             const int n_ener0,
//                             const double *parameter,
//                             int ifl,
//                             double *photar,
//                             double *photer,
//                             const char *init) {
//
//  const int n_parameter = 15;
//  int status = EXIT_SUCCESS;
//  tdrelxilllpion_nthcomp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxilllpionCp model failed", &status);
//}
//
//// XSpec local model function (defined in the lmodel*.dat file)
//void lmodxillverdensnthcomp(const double *ener0,
//                            const int n_ener0,
//                            const double *parameter,
//                            int ifl,
//                            double *photar,
//                            double *photer,
//                            const char *init) {
//
//  const int n_parameter = 8;
//  int status = EXIT_SUCCESS;
//  tdxillverdens_nthcomp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating xillverDCp model failed", &status);
//}
//
//// XSpec local model function (defined in the lmodel*.dat file)
//void lmodrelxilllpdensnthcomp(const double *ener0,
//                              const int n_ener0,
//                              const double *parameter,
//                              int ifl,
//                              double *photar,
//                              double *photer,
//                              const char *init) {
//
//  const int n_parameter = 15;
//  int status = EXIT_SUCCESS;
//  tdrelxilllpdens_nthcomp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxilllpDCp model failed", &status);
//}
//
//// XSpec local model function (defined in the lmodel*.dat file)
//void lmodrelxilldensnthcomp(const double *ener0,
//                            const int n_ener0,
//                            const double *parameter,
//                            int ifl,
//                            double *photar,
//                            double *photer,
//                            const char *init) {
//
//  const int n_parameter = 14;
//  int status = EXIT_SUCCESS;
//  tdrelxilldens_nthcomp(ener0, n_ener0, photar, parameter, n_parameter, &status);
//
//  if (status != EXIT_SUCCESS)
//    RELXILL_ERROR("evaluating relxill_dens model failed", &status);
//}
//
