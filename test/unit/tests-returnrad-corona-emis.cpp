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


// Unit tests of the emissivity profile and line profile of the returning radiation of the corona (lamp post models,
// Relreturn_Corona.cpp). Run with bin/tests "[rrad-corona][emis]".

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "Rellp.h"
#include "Relbase.h"
#include "Relprofile.h"
#include "Relreturn_Table.h"
#include "writeOutfiles.h"

extern "C" {
#include "relutility.h"
#include "envoptions.h"
}

static void invert_emis_profile(emisProfile *emis) {
  invertArray(emis->re, emis->nr);
  invertArray(emis->emis, emis->nr);

  if (emis->del_inc != nullptr) {
    invertArray(emis->del_inc, emis->nr);
  }

  if (emis->del_emit != nullptr) {
    invertArray(emis->del_emit, emis->nr);
  }
}

static void write_emis_profile(const std::string &fname, emisProfile *emis_profile) {
  write_data_to_file(fname.c_str(), emis_profile->re, emis_profile->emis, emis_profile->nr);
}

TEST_CASE (

"corona returning radiation: emissivity profile rebinned from the table grid agrees with the fine grid within 5%"
,
"[rrad-corona][emis]"
)
 {

  int status = EXIT_SUCCESS;

  double precRebinCoarse = 0.05;

  LocalModel lmod(ModelName::relline_lp);
  lmod.set_par(XPar::switch_switch_returnrad, 0);  // need this for the comparison (which is without return_rad)
  const relParam* rel_param = lmod.rel_params();

  RelSysPar *sysPar = get_system_parameters(rel_param, &status);

  returningFractions *dat = get_rrad_fractions(rel_param->a, rel_param->rin, rel_param->rout, &status);
  auto* rad = new double[dat->nrad];
  auto nrad = dat->nrad;
  std::memcpy(rad, dat->rad.data(), nrad * sizeof(double));  // one-line copy
  delete_returningFractions(&dat);

 invertArray(rad, nrad);
 emisProfile *emisCoarse = calc_emis_profile(rad, nrad, rel_param, &status);
 invert_emis_profile(emisCoarse);

  emisProfile *emisFineReference = calc_emis_profile(sysPar->re, sysPar->nr, rel_param, &status);
  emisProfile *emisRebin = new_emisProfile(sysPar->re, sysPar->nr, &status);
  rebin_emisprofile_on_radial_grid(emisRebin, emisCoarse, &status);

  REQUIRE(status==EXIT_SUCCESS);

  write_emis_profile("test_emisCoarse.dat", emisCoarse);
  write_emis_profile("test_emisFineReference.dat", emisFineReference);
  write_emis_profile("test_emisRebin.dat", emisRebin);

  for (int ii = 20; ii < emisRebin->nr;  ii++) {  // last bin in coarse grid deviates, so skip, as it is not the interpolation
    if (is_debug_run()) {
      printf(" rad: %.3e :  %e (ref=%e, ratio=%e) \n",
             emisRebin->re[ii], emisRebin->emis[ii], emisFineReference->emis[ii],
             emisRebin->emis[ii] / emisFineReference->emis[ii]);
    }
    REQUIRE(fabs(emisRebin->emis[ii] / emisFineReference->emis[ii] - 1) < precRebinCoarse);
  }

    free_emisProfile(emisCoarse);
    free_emisProfile(emisFineReference);
    free_emisProfile(emisRebin);
    delete[] rad;
}

TEST_CASE (

"corona returning radiation: relline_lp line profile changes with returning radiation"
,
"[rrad-corona][emis]"
)
 {

  int status = EXIT_SUCCESS;

  LocalModel lmod(ModelName::relline_lp);
  relParam* rel_param = lmod.get_rel_params_struct();
  rel_param->a = 0.998;
  rel_param->height = 5.0;
  rel_param->rout = 1000.0;

  const char* env_name_outfiles = "RELXILL_WRITE_OUTFILES";

  TestSpectrum *spec = getNewSpec(0.05, 10, 1000, &status);

  rel_param->return_rad=0;
  relline_spec_multizone *rel_profile_norrad = relbase(spec->ener, spec->nbins, rel_param, &status);

  setenv(env_name_outfiles, "1", 1);
  rel_param->return_rad = 1;
  relline_spec_multizone *rel_profile = relbase(spec->ener, spec->nbins, rel_param, &status);
  unsetenv(env_name_outfiles);

  REQUIRE(status==EXIT_SUCCESS);
  REQUIRE(rel_profile->n_zones == 1);
  REQUIRE(calcSum(rel_profile->flux[0], rel_profile->n_ener) > 1e-8);

  double abs_diff_profiles[rel_profile->n_ener];
  for (int ii = 0; ii < rel_profile->n_ener; ii++) {
    abs_diff_profiles[ii] = abs(rel_profile_norrad->flux[0][ii] - rel_profile->flux[0][ii]);
  }

  INFO(" require that the line profile with and without return radiation differs ");
  CHECK_THAT(calcSum(abs_diff_profiles, rel_profile->n_ener),
             !Catch::Matchers::WithinAbs(0.0, 1e-2)
  );

  free_Spectrum(&spec);
  delete_rel_params(rel_param);
}

TEST_CASE (

"corona returning radiation: emissivity profile is calculated (written to __output_emis_profile_rrad.dat)"
,
"[rrad-corona][emis]"
)
 {

  int status = EXIT_SUCCESS;

  LocalModel lmod(ModelName::relline_lp);
  lmod.set_par(XPar::switch_switch_returnrad, -1);  // need this for the comparison (which is without return_rad)
  const relParam *rel_param = lmod.rel_params();
  RelSysPar *sysPar = get_system_parameters(rel_param, &status);
  emisProfile *emis_profile = calc_emis_profile(sysPar->re, sysPar->nr, rel_param, &status);
  write_emis_profile("__output_emis_profile_rrad.dat", emis_profile);

  REQUIRE(status == EXIT_SUCCESS);

  free_emisProfile(emis_profile);
}
