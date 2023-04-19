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
#include "IonGradient.h"
#include "common-functions.h"

#define PREC 1e-6


TEST_CASE(" Ion Grad PL Index ", "[iongrad]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  int status = EXIT_SUCCESS;

  LocalModel local_model{ModelName::relxilllpCp};

  local_model.eval_model(spec);

  xillParam *xill_param = local_model.get_xill_params();
  relParam *rel_param = local_model.get_rel_params();

  // relline_spec_multizone *rel_profile = relbase(spec.energy, spec.num_flux_bins(), rel_param, &status);
  RelSysPar *sys_par = get_system_parameters(rel_param, &status);

  RadialGrid radial_grid{rel_param->rin, rel_param->rout, rel_param->num_zones, rel_param->height};
  IonGradient ion_gradient{radial_grid, rel_param->ion_grad_type, xill_param->iongrad_index};
  ion_gradient.calculate_gradient(*(sys_par->emis), PrimarySourceParameters{local_model.get_model_params()});

  for (int ii = 0; ii < rel_param->num_zones; ii++) {
    REQUIRE(ion_gradient.dens[ii] >= 15.0);
    REQUIRE(ion_gradient.dens[ii] <= 22.0);
    REQUIRE(ion_gradient.lxi[ii] >= 0.0);
    REQUIRE(ion_gradient.lxi[ii] <= 4.7);
  }

}


TEST_CASE(" Test Alpha Model (writing output) ", "[iongrad-write]") {

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  LocalModel local_model{ModelName::relxilllpCp};

  const char *env_outfiles = "RELXILL_WRITE_FILES";
  setenv(env_outfiles, "1", 1);
  setenv("RELXILL_NUM_RZONES", "50", 1);

  local_model.set_par(XPar::logn, 15.0);
  local_model.set_par(XPar::h, 3);
  local_model.set_par(XPar::a, 0.99);
  local_model.set_par(XPar::switch_iongrad_type, 2);
  local_model.eval_model(spec);

  local_model.set_par(XPar::a, 0.998);
  local_model.eval_model(spec);

  unsetenv(env_outfiles);
}

TEST_CASE(" Exec single iongrad model, should change with xindex", "[iongrad]") {
  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxilllpCp);

  auto spec = default_spec.get_xspec_spectrum();
  lmod.set_par(XPar::switch_iongrad_type, 1);
  lmod.set_par(XPar::iongrad_index, 1.0);

  REQUIRE_NOTHROW(lmod.eval_model(spec));
  double sum = sum_flux(spec.flux, spec.num_flux_bins());

  REQUIRE_NOTHROW(lmod.eval_model(spec));

  auto spec2 = default_spec.get_xspec_spectrum();
  lmod.set_par(XPar::iongrad_index, 2.0);

  REQUIRE_NOTHROW(lmod.eval_model(spec2));
  double sum2 = sum_flux(spec2.flux, spec2.num_flux_bins());

  REQUIRE(sum > 1e-8);
  REQUIRE(sum2 > 1e-8);

  REQUIRE( fabs(sum - sum2) > 1e-8);

}


IonGradient get_ion_gradient(LocalModel const &lmod){
  xillParam *xill_param = lmod.get_xill_params();
  relParam *rel_param = lmod.get_rel_params();

  // relline_spec_multizone *rel_profile = relbase(spec.energy, spec.num_flux_bins(), rel_param, &status);
  int status = EXIT_SUCCESS;
  RelSysPar *sys_par = get_system_parameters(rel_param, &status);

  RadialGrid const radial_grid{rel_param->rin, rel_param->rout, rel_param->num_zones, rel_param->height};

  IonGradient ion_gradient{radial_grid, rel_param->ion_grad_type, xill_param->iongrad_index};
  ion_gradient.calculate_gradient(*(sys_par->emis), PrimarySourceParameters{lmod.get_model_params()});

  delete rel_param;
  delete xill_param;

  return ion_gradient;
}

TEST_CASE(" Test Constant Density by Env Variable ", "[iongrad-const-density]") {

  DefaultSpec const default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel lmod{ModelName::relxilllpAlpha};

  const char* env_density = "RELXILL_CONSTANT_DENSITY";
  setenv(env_density, "1", 1);
  setenv("RELXILL_NUM_RZONES","10",1);

  const double dens_param = 16;
  lmod.set_par(XPar::logn, dens_param);
  lmod.eval_model(spec);
  auto ion_grad = get_ion_gradient(lmod);
  const int num_zones = ion_grad.radial_grid.num_zones();

  unsetenv(env_density);

  lmod.set_par(XPar::logn, dens_param*1.1);
  lmod.eval_model(spec);
  lmod.set_par(XPar::logn, dens_param*1.2);
  lmod.eval_model(spec);


  for (int ii = 0; ii < num_zones; ii++) {
    REQUIRE( fabs(dens_param - ion_grad.dens[ii]) < 1e-6 );
  }
}

TEST_CASE(" Test Density Gradient ", "[iongrad-density]") {

  DefaultSpec const default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel lmod{ModelName::relxilllpAlpha};

  setenv("RELXILL_NUM_RZONES", "50", 1);

  const double dens_param = 16;
  lmod.set_par(XPar::logn, dens_param);
  // lmod.eval_model(spec);
  auto ion_grad = get_ion_gradient(lmod);
  const int num_zones = ion_grad.radial_grid.num_zones();

  // require that the minimal density is given by the dens_param
  int index_dens_min = 0;
  double dens_min = ion_grad.dens[index_dens_min];
  for (int ii = 0; ii < ion_grad.nzones(); ii++) {
    if (ion_grad.dens[ii] < dens_min) {
      dens_min = ion_grad.dens[ii];
      index_dens_min = ii;
    }
  }

  // value where the density is minimal
  const double r_dens_min = (25. / 9.) * ion_grad.radial_grid.radius[0];  // 25/9*rin

  // minimal density at the correct radius?
  REQUIRE(ion_grad.radial_grid.radius[index_dens_min] < r_dens_min);
  REQUIRE(ion_grad.radial_grid.radius[index_dens_min + 1] > r_dens_min);

  // value of the minimal density as expected?
  REQUIRE(fabs(dens_param - dens_min) < 1e-2);
  REQUIRE(fabs(dens_param - ion_grad.dens[0]) > 1e-2);  // for the given parameters the inner zone should be larger
}

TEST_CASE(" Test Density Gradient ", "[iongrad-valgrind]") {

  DefaultSpec const default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  LocalModel lmod{ModelName::relxilllpAlpha};

  setenv("RELXILL_NUM_RZONES","20",1);

  const double dens_param = 16;
  const double dens_param_delta = 1;
  lmod.set_par(XPar::logn, dens_param);
  lmod.eval_model(spec);

  const char* env_density = "RELXILL_CONSTANT_DENSITY";
  setenv(env_density, "1", 1);

  const int num_eval = 10;
  for (int ii = 1; ii < num_eval; ii++) {
    lmod.set_par(XPar::logn, dens_param+ dens_param_delta*float(ii)/float(num_eval));
    lmod.eval_model(spec);
  }

  unsetenv(env_density);
}
