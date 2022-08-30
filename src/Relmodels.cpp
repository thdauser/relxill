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

#include "Relmodels.h"

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
    printf(" *** relxill error : Spin a > 0.9982, model evaluation failed (value is %f) \n", param->a);
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
        warned_height = 1;
      }
      param->height = r_event * h_fac;
    }
  }

}


/** BASIC XILLVER MODEL FUNCTION **/
void xillver_base(double *ener_inp, const int n_ener0, double *photar, const ModelParams &inp_params, int *status) {

  xillParam *xill_param = get_xill_params(inp_params);

  // call the function which calculates the xillver spectrum
  xillSpec *spec = get_xillver_spectra(xill_param, status);
  CHECK_STATUS_VOID(*status);

  // =4= rebin to the input grid
  assert(spec->n_incl == 1); // make sure there is only one spectrum given (for the chosen inclination)

  /** add the dependence on incl, assuming a semi-infinite slab **/
  norm_xillver_spec(spec, xill_param->incl);

  rebin_spectrum(ener_inp, photar, n_ener0, spec->ener, spec->flu[0], spec->n_ener);

  add_primary_component(ener_inp, n_ener0, photar, nullptr, xill_param, nullptr, status);

  free_xill_spec(spec);
  delete xill_param;
}

void relline_base(double *ener1keV, double *photar, const int n_ener, relParam *param_struct, int *status) {

  relline_spec_multizone *spec = relbase(ener1keV, n_ener, param_struct, status);

  for (int ii = 0; ii < n_ener; ii++) {
    photar[ii] = spec->flux[0][ii];
  }
}

