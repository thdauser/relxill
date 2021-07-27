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

  rel_spec *rel_profile = get_stdRelProfile(&status);

  const double ReferenceRelatStdFlux = 8.371512e-01;
  double relatFlux = calc_RelatFluxInStdBand(rel_profile);
  compareReferenceFlux(relatFlux, ReferenceRelatStdFlux, &status);
  REQUIRE( status == EXIT_SUCCESS );

}

TEST_CASE(" compare standard xillver evaluation with reference flux") {

  int status = EXIT_SUCCESS;

  rel_spec *rel_profile = nullptr;
  double *xill_spec = nullptr;
  init_std_relXill_spec(&rel_profile, &xill_spec, &status);

  const double ReferenceXillverStdFlux = 1.802954e+01;
  double xillFlux = calc_XillverFluxInStdBand(xill_spec, rel_profile->ener, rel_profile->n_ener);
  compareReferenceFlux(xillFlux, ReferenceXillverStdFlux, &status);
  REQUIRE( status == EXIT_SUCCESS );

}

