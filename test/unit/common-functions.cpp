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

#include "common-functions.h"

double sum_flux(const double *flux, int nbins) {

  double sum = 0.0;
  for (int ii = 0; ii < nbins; ii++) {
    sum += flux[ii];
  }
  return sum;
}

void eval_xspec_lmod_default(ModelName model_name, const DefaultSpec& default_spec){

  const double *xspec_parameters = get_xspec_default_parameter_array(model_name);

  xspec_C_wrapper_eval_model(model_name,
                               xspec_parameters,
                               default_spec.flux,
                               default_spec.num_flux_bins,
                               default_spec.energy);

  delete[] xspec_parameters;

}



void get_RelProfileConstEmisZones(relline_spec_multizone **p_rel_profile, relParam **p_rel_param, int nzones, int *status) {

  // relline is per default re-normalized, but we need to test the physical normalization here
  const char* env_physnorm = "RELLINE_PHYSICAL_NORM";
  setenv(env_physnorm, "1", 1);

  LocalModel lmod(ModelName::relline);
  lmod.set_par(XPar::index1, 0.0);
  lmod.set_par(XPar::index2, 0.0);
  lmod.set_par(XPar::a, 0.0);
  lmod.set_par(XPar::rout, 1000);
 // lmod.set_par(XPar::limb, 1.0);

  int n_ener;
  double *ener;
  get_relxill_conv_energy_grid(&n_ener, &ener, status);

  *p_rel_param = lmod.get_rel_params();
  (*p_rel_param)->num_zones = nzones;

  assert(*p_rel_param != nullptr);
  *p_rel_profile = relbase(ener, n_ener, *p_rel_param, status);
  setenv(env_physnorm, "0", 1);

}

relline_spec_multizone *get_stdRelProfile(int *status) {

  LocalModel lmod(ModelName::relline);

  int n_ener;
  double *ener;
  get_relxill_conv_energy_grid(&n_ener, &ener, status);

  return relbase(ener, n_ener, lmod.get_rel_params(), status);
}

xillSpec *get_std_xill_spec(int *status) {
  LocalModel lmod(ModelName::xillver);
  xillSpec *xill_spec = get_xillver_spectra(lmod.get_xill_params(), status);
  return xill_spec;
}


/**
 * calculate a relativistic profile and a xillver spectrum for standard parameters on the
 * same energy grid (for testing)
 * @param rel_profile [output]
 * @param xill_spec_output [output]
 * @param status
 */
void init_std_relXill_spec(relline_spec_multizone **rel_profile, double **xill_spec_output, int *status) {

  CHECK_STATUS_VOID(*status);

  *rel_profile = get_stdRelProfile(status);

  auto *xill_flux = new double[(*rel_profile)->n_ener];
  CHECK_MALLOC_VOID_STATUS(xill_flux, status)

  xillSpec *xill_spec_table = get_std_xill_spec(status);
  rebin_spectrum((*rel_profile)->ener,
                 xill_flux,
                 (*rel_profile)->n_ener,
                 xill_spec_table->ener,
                 xill_spec_table->flu[0],
                 xill_spec_table->n_ener);
  *xill_spec_output = xill_flux;

}

double calc_FluxInStdBand(const double *flux, double *ener, int n) {

  const double ELO_REF_RELPROFILE = 0.1;
  const double EHI_REF_RELPROFILE = 1.0;

  return calcSumInEnergyBand(flux, n, ener, ELO_REF_RELPROFILE, EHI_REF_RELPROFILE);
}

double calc_RelatFluxInStdBand(const relline_spec_multizone *spec) {

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
