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


// Unit tests of the correction factors of the returning radiation of the corona (flux and energy-shift corrections
// from the xillver spectra, Relreturn_Corona.cpp and Relxill.cpp). Run with bin/tests "[rrad-corona][corrfactors]".

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "Xillspec.h"
#include "Xilltable_handling.h"
#include "Relxill.h"
#include "Relreturn_Corona.h"
#include "Relreturn_Table.h"
#include "writeOutfiles.h"

extern "C" {
#include "relutility.h"
#include "xilltable_reading.h"
}

TEST_CASE (

"correction factors: flux correction is in (0.5, 1) and larger for high than for low ionization"
,
"[rrad-corona][corrfactors]"
)
 {
  //

  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  lmod.set_par(XPar::logxi, 3.0);
  double flux_corr1=0;
  xillSpec* xill_spec = get_xillver_spectra(lmod.xill_params(), &status);
    auto xilltab_param = get_xilltab_param(lmod.xill_params(), &status);
  get_xillver_fluxcorrection_factors(xill_spec, &flux_corr1, nullptr,
                                     xilltab_param,
                                     &status);
  REQUIRE(flux_corr1 > 0.5);
  REQUIRE(flux_corr1 < 1.0); // only works for logxi<=4

  lmod.set_par(XPar::logxi, 0.0);
  xillSpec* xill_spec2 = get_xillver_spectra(lmod.xill_params(), &status);
  double flux_corr2=0;
    auto xilltab_param2 = get_xilltab_param(lmod.xill_params(), &status);
  get_xillver_fluxcorrection_factors(xill_spec2, &flux_corr2, nullptr,
                                     xilltab_param2 ,
                                     &status);
  free(xilltab_param);
  free(xilltab_param2);

  REQUIRE(flux_corr1 > 2 * flux_corr2);  // for low instead high ionization it is much lower, only works if logxi_hi>3.5

  free_xill_spec(xill_spec);
  free_xill_spec(xill_spec2);

}

TEST_CASE (

"correction factors: gshift correction is above 1 for logxi=3 and below 1 for logxi=0"
,
"[rrad-corona][corrfactors]"
)
 {
  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  double gshift_ref_value = 1.5; // correction factor calculated for this value
  double gam = 2.0;

  lmod.set_par(XPar::gamma, gam);
  lmod.set_par(XPar::logxi, 3.0);
  double gshift_flux_corr1;
  xillSpec* xill_spec = get_xillver_spectra(lmod.xill_params(), &status);
    auto xilltab_param = get_xilltab_param(lmod.xill_params(), &status);
  get_xillver_fluxcorrection_factors(xill_spec, nullptr, &gshift_flux_corr1, xilltab_param , &status);
    free(xilltab_param);

  REQUIRE(gshift_flux_corr1 > 1);
  REQUIRE(gshift_flux_corr1 < pow(gshift_ref_value,gam+0.5)/pow(gshift_ref_value,gam));

  lmod.set_par(XPar::logxi, 0.0);
  double gshift_flux_corr2;
  xillSpec* xill_spec2 = get_xillver_spectra(lmod.xill_params(), &status);
    auto xilltab_param2 = get_xilltab_param(lmod.xill_params(), &status);
  get_xillver_fluxcorrection_factors(xill_spec2, nullptr, &gshift_flux_corr2,
                                     xilltab_param2, &status);
    free(xilltab_param2);

  REQUIRE(gshift_flux_corr2 < 1 );
  REQUIRE(gshift_flux_corr1 > 1.5 * gshift_flux_corr2);  // for low instead high ionization it is much lower, only works if logxi_hi>3.5

  free_xill_spec(xill_spec);
  free_xill_spec(xill_spec2);

}

TEST_CASE (

"correction factors: corrected_gshift_fluxboost_factor is above (below) g^gamma for a factor above (below) 1"
,
"[rrad-corona][corrfactors]"
)
 {

  double gamma = 2.0;

  double xill_gshift_fac = 1.2;
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) > pow(1.5, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) > pow(0.3, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) > 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) < 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) > 0);

  xill_gshift_fac = 0.9;
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) < pow(1.5, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) < pow(0.3, gamma));
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 1.5, gamma) > 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) < 1);
  REQUIRE(corrected_gshift_fluxboost_factor(xill_gshift_fac, 0.3, gamma) > 0);

}

TEST_CASE (

"correction factors: flux and gshift correction vs. logxi (written to test-*-corr-returnrad_gam*.dat)"
,
"[rrad-corona][corrfactors]"
)
 {
  //

  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  const int n_logxi = 100;
  double logxi[n_logxi] = {0};
  double flux_corr[n_logxi] = {0};
  double gshift_corr[n_logxi] = {0};

  double logxi_min = 0.0;
  double logxi_max = 4.0;

  const int n_gam = 3;
  double values_gamma[n_gam] = {1.5, 2.0, 2.5};

  char buffer[100];
  for (double gam : values_gamma) {
    lmod.set_par(XPar::gamma, gam);
    for (int ii = 0; ii < n_logxi; ii++) {
      logxi[ii] = (logxi_max - logxi_min) * (ii * 1.0 / (n_logxi - 1)) + logxi_min;
      lmod.set_par(XPar::logxi, logxi[ii]);
      xillSpec* xill_spec = get_xillver_spectra(lmod.xill_params(), &status);
        auto xilltab_param = get_xilltab_param(lmod.xill_params(), &status);
      get_xillver_fluxcorrection_factors(xill_spec, &(flux_corr[ii]), &(gshift_corr[ii]),
                                         xilltab_param, &status);
      free(xilltab_param);
      free_xill_spec(xill_spec);
    }
    int retval = sprintf(buffer, "test-flux-corr-returnrad_gam%.1f.dat",gam);
    write_data_to_file(buffer, logxi, flux_corr, n_logxi);
    retval = sprintf(buffer, "test-gshift-corr-returnrad_gam%.1f.dat",gam);
    write_data_to_file(buffer, logxi, gshift_corr, n_logxi);
  }

}

static void require_all_values_the_same(const double *val, const int n) {
  double avg = calcSum(val, n) / n;
  REQUIRE(fabs(avg - val[0]) < 1e-6);
}

static void require_values_differ(const double *val, const int n) {
  double avg = calcSum(val, n) / n;
  REQUIRE(fabs((avg - val[0]) / avg) > 1e-2);
  REQUIRE(fabs((avg - val[n - 1]) / avg) > 1e-2);
}

TEST_CASE (

"correction factors: constant over the disk for a constant spectrum, varying for an ionization gradient"
,
"[rrad-corona][corrfactors]"
)
 {

  int status = EXIT_SUCCESS;

  LocalModel lmod{ModelName::relxilllp};
  lmod.set_par(XPar::switch_switch_returnrad, 1);

  relParam* rel_param = lmod.get_rel_params_struct();

  const auto n_zones = rel_param->num_zones;
  double logxi_min = 0.0;
  double logxi_max = 4.0;

  auto radial_grid = RadialGrid(rel_param->rin, rel_param->rout, n_zones, rel_param->height);

  xillTableParam* xill_table_param[n_zones];
  xillSpec* xill_spec[n_zones];

  // TEST contant spectrum over the disk
  for (int ii = 0; ii < n_zones; ii++) {
    xill_table_param[ii] = get_xilltab_param(lmod.xill_params(), &status);
    xill_spec[ii] = get_xillver_spectra_table(xill_table_param[ii], &status);
  }

  rel_param->rrad_corr_factors =
      calc_rrad_corr_factors(xill_spec, radial_grid, xill_table_param, &status);

    require_all_values_the_same(rel_param->rrad_corr_factors->corrfac_flux, n_zones);
  require_all_values_the_same(rel_param->rrad_corr_factors->corrfac_gshift, n_zones);

  returningFractions *ret_fractions = get_rrad_fractions(rel_param->a, rel_param->rin, rel_param->rout, &status);

  rradCorrFactors* tab_corrfactors =
      rebin_corrfactors_to_rradtable_grid(rel_param->rrad_corr_factors, ret_fractions, &status);

  require_all_values_the_same(tab_corrfactors->corrfac_flux, tab_corrfactors->n_zones);
  require_all_values_the_same(tab_corrfactors->corrfac_gshift, tab_corrfactors->n_zones);

  // as all values are the same, re-binning should not change the actual value
  REQUIRE(rel_param->rrad_corr_factors->corrfac_flux[0] == tab_corrfactors->corrfac_flux[0]);
  REQUIRE(rel_param->rrad_corr_factors->corrfac_gshift[0] == tab_corrfactors->corrfac_gshift[0]);

  free_rrad_corr_factors(&tab_corrfactors);

  // TEST ionization gradient
  xillSpec* xill_spec_grad[n_zones];
  for (int ii = 0; ii < n_zones; ii++) {
    xill_table_param[ii]->lxi = (logxi_max - logxi_min) * (ii * 1.0 / (n_zones - 1)) + logxi_min;
    xill_spec_grad[ii] = get_xillver_spectra_table(xill_table_param[ii], &status);
  }

  rradCorrFactors* grad_corrfactors =
      calc_rrad_corr_factors(xill_spec_grad, radial_grid, xill_table_param, &status);

  require_values_differ(grad_corrfactors->corrfac_flux, n_zones);
  require_values_differ(grad_corrfactors->corrfac_gshift, n_zones);

  rradCorrFactors* tab_grad_corrfactors =
      rebin_corrfactors_to_rradtable_grid(grad_corrfactors, ret_fractions, &status);
  require_values_differ(tab_grad_corrfactors->corrfac_flux, tab_grad_corrfactors->n_zones);
  require_values_differ(tab_grad_corrfactors->corrfac_gshift, tab_grad_corrfactors->n_zones);

    free_rrad_corr_factors(&tab_grad_corrfactors);
    free_rrad_corr_factors(&grad_corrfactors);
    for (int ii = 0; ii < n_zones; ii++)
    {
        free_xill_spec(xill_spec[ii]);
        free_xill_spec(xill_spec_grad[ii]);
        free(xill_table_param[ii]);
        xill_table_param[ii] = nullptr;
    }
    delete_rel_params(rel_param);
    delete_returningFractions(&ret_fractions);
}
