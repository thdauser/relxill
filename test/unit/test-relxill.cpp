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
#include "common-functions.h"

extern "C" {
#include "writeOutfiles.h"
}



TEST_CASE(" compare standard relline profile with reference flux "){

  int status = EXIT_SUCCESS;

  relline_spec_multizone *rel_profile = get_stdRelProfile(&status);

  const double ReferenceRelatStdFlux = 8.533758e-01;
  double relatFlux = calc_RelatFluxInStdBand(rel_profile);
  compareReferenceFlux(relatFlux, ReferenceRelatStdFlux, &status);
  REQUIRE( status == EXIT_SUCCESS );

}

TEST_CASE(" test the relxill normalization", "[relxill-norm]") {

  auto default_spec = DefaultSpec(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, 3000);
  auto spec = default_spec.get_xspec_spectrum();


  double ref_value = 3.0; // normalized at 1keV to have 1 cts/sec/keV/cm2 at 1keV
  int ind = binary_search(spec.energy, spec.num_flux_bins(), ref_value);
  double dE_ref = spec.energy[ind+1] - spec.energy[ind];


  LocalModel lmod(ModelName::relxilllp);
  lmod.set_par(XPar::h, 30);
  lmod.eval_model(spec);

  // by default relxill should not be normalized
  REQUIRE(  fabs(spec.flux[ind]/dE_ref - 1 ) > 1e-2);

  // now it is normalized
  setenv("RELXILL_RENORMALIZE","1",1);
  lmod.eval_model(spec);
  REQUIRE(  fabs(spec.flux[ind]/dE_ref -1 ) < 1e-2);

  // should stay the same even if we change a parameter like the height
  lmod.set_par(XPar::h, 20);
  lmod.eval_model(spec);
  REQUIRE(  fabs(spec.flux[ind]/dE_ref - 1) < 1e-2);


  setenv("RELXILL_RENORMALIZE","0",1);
}

