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

#include "catch2/catch_amalgamated.hpp"

#include "LocalModel.h"
#include "xspec_wrapper_lmodels.h"
#include "XspecSpectrum.h"
#include "common-functions.h"
#include "Rellp.h"

#define PREC 1e-6

TEST_CASE(" test printing of reflection strength", "[rellp]" ){

  DefaultSpec default_spec{};
  LocalModel lmod(ModelName::relxilllp);

  const char *env = "RELXILL_PRINT_DETAILS";

  lmod.set_par(XPar::refl_frac, 1);
  lmod.set_par(XPar::switch_switch_reflfrac_boost, 1);
  auto spec = default_spec.get_xspec_spectrum();

  lmod.eval_model(spec);

  setenv(env, "1", 1);

  lmod.eval_model(spec);

  unsetenv(env);

}

TEST_CASE(" beta>0 of a lamp post changes the primary spectrum", "[beta]") {

  DefaultSpec default_spec{};
  LocalModel lmod(ModelName::relxilllpCp);

  lmod.set_par(XPar::refl_frac, 0);
  lmod.set_par(XPar::switch_switch_reflfrac_boost, 1);
  lmod.set_par(XPar::kte, 100);
  lmod.set_par(XPar::h, 30);
  lmod.set_par(XPar::incl, 30);
  auto spec = default_spec.get_xspec_spectrum();

  lmod.set_par(XPar::beta, 0.66);
  lmod.eval_model(spec);

  double sum0 = calcSumInEnergyBand(spec.flux, spec.num_flux_bins(), spec.energy, 1, 10);

  lmod.set_par(XPar::beta, 0.0);
  lmod.eval_model(spec);

  double sum1 = calcSumInEnergyBand(spec.flux, spec.num_flux_bins(), spec.energy, 1, 10);

  REQUIRE(sum0 > sum1);

}

TEST_CASE(" normalization of LP primary spectrum", "[beta]") {

  int status = EXIT_SUCCESS;

  EnerGrid *egrid = get_coarse_xillver_energrid(&status);

  DefaultSpec default_spec{};
  LocalModel lmod(ModelName::relxilllpCp);
  lmod.set_par(XPar::refl_frac, 0);
  lmod.set_par(XPar::switch_switch_reflfrac_boost, 1);
  lmod.set_par(XPar::kte, 100);
  lmod.set_par(XPar::h, 30);
  lmod.set_par(XPar::incl, 20);
  xillTableParam *xill_param = get_xilltab_param(lmod.get_xill_params(), &status);

  auto prim_spec_source = new double[egrid->nbins];
  xill_param->ect = 100;
  calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins, xill_param, &status);
  double norm_fac1 = 1./calcNormWrtXillverTableSpec(prim_spec_source, egrid->ener, egrid->nbins, &status);

  xill_param->ect = 100;
  calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins, xill_param, &status, 0.3);
  double norm_fac2 = 1. / calcNormWrtXillverTableSpec(prim_spec_source, egrid->ener, egrid->nbins, &status);

  REQUIRE(norm_fac1 < norm_fac2);

  // printf(" norm1=%e   norm2=%e\n", norm_fac1, norm_fac2);

}

TEST_CASE(" Change of Ecut on the disk with beta>0  ", "[beta]") {

  int status = EXIT_SUCCESS;

  double kte_primary = 40.0;

  LocalModel local_model{ModelName::relxilllpCp};
  local_model.set_par(XPar::switch_iongrad_type, 0);
  local_model.set_par(XPar::h, 50.0);
  local_model.set_par(XPar::kte, kte_primary);

  xillParam *xill_param = local_model.get_xill_params();
  relParam *rel_param = local_model.get_rel_params();
  RelSysPar *sys_par = get_system_parameters(rel_param, &status);

  RadialGrid radial_grid{rel_param->rin, rel_param->rout, rel_param->num_zones, rel_param->height};
  IonGradient ion_gradient{radial_grid, rel_param->ion_grad_type};
  ion_gradient.calculate_gradient(*(sys_par->emis), rel_param, xill_param);

  double ecut_in_0 = ion_gradient.get_ecut_disk_zone(rel_param, kte_primary, 0);
  double ecut_out_0 = ion_gradient.get_ecut_disk_zone(rel_param, kte_primary, 1);


  /*  for (int ii=0; ii<rel_param->num_zones; ii++){
      double rad = 0.5 * (ion_gradient.radial_grid.radius[ii] + ion_gradient.radial_grid.radius[ii + 1]);
      printf(" del_emit=%.1f, rad=%.2e, kTe= %.2e \n",
             ion_gradient.del_emit[ii]*180.0/M_PI, rad,
             ion_gradient.get_ecut_disk_zone(rel_param, kte_primary, ii) );
    } */

  // now set the velocity to beta>0
  rel_param->beta = 0.66;
  sys_par = get_system_parameters(rel_param, &status);
  ion_gradient.calculate_gradient(*(sys_par->emis), rel_param, xill_param);
  double ecut_in_beta = ion_gradient.get_ecut_disk_zone(rel_param, kte_primary, 0);



  /* for (int ii=0; ii<rel_param->num_zones; ii++){
     double rad = 0.5 * (ion_gradient.radial_grid.radius[ii] + ion_gradient.radial_grid.radius[ii + 1]);
     printf(" del_emit=%.1f, rad=%.2e, kTe= %.2e \n",
            ion_gradient.del_emit[ii]*180.0/M_PI, rad,
            ion_gradient.get_ecut_disk_zone(rel_param, kte_primary, ii) );
   } */


  // require that the first zone is most blue-shifted
  REQUIRE(ecut_in_0 > ecut_out_0);

  // require that if we have beta>0 ecut is redshift wrt to beta=0
  REQUIRE(ecut_in_0 > ecut_in_beta);

}


/*

static void test_pointSourceDecision(int *status) {

  printf(" *** testing: pointSourcedDecision\n");

  relParam *rel_param = get_std_relParam_relxilllpDCp(status);

  rel_param->height = 3.0;
  rel_param->htop = 3.0;
  assert(modelLampPostPointsource(rel_param) == 1);

  rel_param->htop = 10.0;
  assert(modelLampPostPointsource(rel_param) == 0);

  rel_param->htop = 1.0;
  assert(modelLampPostPointsource(rel_param) == 1);

  rel_param->htop = 0.0;
  assert(modelLampPostPointsource(rel_param) == 1);

  CHECK_RELXILL_DEFAULT_ERROR(status);

}

static void test_extendedGeometryHeight(int *status) {

  printf(" *** testing: extendedGeometryHeight\n");

  relParam *rel_param = get_std_relParam_relxilllpDCp(status);

  rel_param->height = 3.0;
  rel_param->htop = 10.0;

  extPrimSource *source = getExtendedJetGeom(rel_param, status);
  assert(*status == EXIT_SUCCESS);

  // test that hbase and htop are set correctly
  assert(fabs(source->heightArr[0] - rel_param->height) < PREC);
  assert(fabs(source->heightArr[source->nh] - rel_param->htop) < PREC);  // we have nh+1 elements

  // test that the height values are ascending
  for (int ii = 0; ii < source->nh; ii++) {
    assert(source->heightArr[ii] < source->heightArr[ii + 1]);
  }

  // test that values are sensible
  assert(source->heightArr[0] >= 1.0);
  assert(source->heightArr[source->nh] <= 1000.0);

}

static void test_calcEmisProfileLp(int *status) {

  printf(" *** testing: calcEmisProfileLp \n");

  relParam *rel_param = get_std_relParam_relxilllpDCp(status);

  rel_param->height = 3.0;
  rel_param->htop = 3.0;
  RelSysPar *sysPar = get_system_parameters(rel_param, status);
  assert(sysPar->emis->normFactorPrimSpec > 0);

  double sumEmisPoint = sumArray(sysPar->emis->emis, sysPar->emis->nr);
  assert(sumEmisPoint > 0);

  rel_param->height = 3.0;
  rel_param->htop = 10.0;
  assert(modelLampPostPointsource(rel_param) == 0);
  RelSysPar *sysPar2 = get_system_parameters(rel_param, status);

  assert(sysPar2->emis->normFactorPrimSpec > 0);
  double sumEmisExt = sumArray(sysPar2->emis->emis, sysPar->emis->nr);
  assert(sumEmisExt > 0);

  assert(sumEmisExt != sumEmisPoint);

  assert(sysPar2->emis->normFactorPrimSpec != sysPar->emis->normFactorPrimSpec);
  assert(sysPar2->emis->photon_fate_fractions->refl_frac != sysPar->emis->photon_fate_fractions->refl_frac);

  CHECK_RELXILL_DEFAULT_ERROR(status);

}
*/

