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
        param->height = r_event * h_fac;
        warned_height = 1;
      }

    }
  }

}


/** BASIC XILLVER MODEL FUNCTION **/
void xillver_base(double *ener_inp, const int n_ener0, double *photar, xillParam *param_struct, int *status) {


  // call the function which calculates the xillver spectrum
  xillSpec *spec = get_xillver_spectra(param_struct, status);
  CHECK_STATUS_VOID(*status);

  // =4= rebin to the input grid
  assert(spec->n_incl == 1); // make sure there is only one spectrum given (for the chosen inclination)

  /** add the dependence on incl, assuming a semi-infinite slab **/
  norm_xillver_spec(spec, param_struct->incl);

  rebin_spectrum(ener_inp, photar, n_ener0, spec->ener, spec->flu[0], spec->n_ener);

  add_primary_component(ener_inp, n_ener0, photar, nullptr, param_struct, nullptr, status);

//  rebin_spectrum(ener, photar, n_ener, ener, flux, n_ener);

  free_xill_spec(spec);
}

void relline_base(double *ener1keV, double *photar, const int n_ener, relParam *param_struct, int *status) {

  relline_spec_multizone *spec = relbase(ener1keV, n_ener, param_struct, status);

  for (int ii = 0; ii < n_ener; ii++) {
    photar[ii] = spec->flux[0][ii];
  }
}

/* get a new relbase parameter structure and initialize it */
relParam *new_relParam(int model_type, int emis_type, int *status) {
  auto *param = (relParam *) malloc(sizeof(relParam));
  if (param == nullptr) {
    RELXILL_ERROR("memory allocation failed", status);
    return nullptr;
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

/* get a new relbase parameter structure and initialize it */
xillParam *new_xillParam(int model_type, int prim_type, int *status) {
  xillParam *param = (xillParam *) malloc(sizeof(xillParam));
  if (param == nullptr) {
    RELXILL_ERROR("memory allocation failed", status);
    return nullptr;
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
  param->boost = -1;
  param->dens = 15;  // the standard value for every table, given in "log-units"
  param->frac_pl_bb = PARAM_DEFAULT;
  param->kTbb = PARAM_DEFAULT;
  param->ion_grad_type = ION_GRAD_TYPE_CONST; // no ion grad
  param->iongrad_index = PARAM_DEFAULT;

  return param;
}
