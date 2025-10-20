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

    Copyright 2025 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "XspecSpectrum.h"
#include "common-functions.h"
#include "Xillspec.h"
#include "Relphysics.h"
#include "Relreturn_Datastruct.h"
#include "Relreturn_Table.h"
#include "Relxill_errorhandling.h"

extern "C" {
#include "relutility.h"
#include "writeOutfiles.h"
#include "xilltable_reading.h"
}


#define PREC 1e-6

TEST_CASE(" Test Loading Single Xillver from RTR Table should fail ", "[rrad-bbody]")
{
    // by default, xill_param->radzone == -1 should fail loading the table
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.xill_params();

    REQUIRE_THROWS(get_xillver_spectra(xill_param, &status));
}


// ------- //
TEST_CASE(" Test Loading Xillver RTR Table ", "[rrad-bbody-vg]")
{
    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.xill_params();
    auto rel_param = lmod.rel_params();

    int status = EXIT_SUCCESS;
    xillTable* xill_tab = nullptr;
    REQUIRE_NOTHROW([&]()
        {
        auto tmp = get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, &status);
        }()
    );

    xill_param->radzone = 0;
    xillSpec* xillspec = nullptr;
    xillspec = get_xillver_spectra(xill_param, &status);
    REQUIRE(xillspec != nullptr);

    // for (int ii = 0; ii < xillspec->n_incl; ii++)
    int ii = 0;
    {
        INFO("Failing index: " << ii);
        double sum = sum_flux(xillspec->flu[ii], xillspec->n_ener);
        printf(" sum: %.5e\n", sum);
        // REQUIRE(sum_flux(xillspec->flu[ii], xillspec->n_ener) > 1e-10);
    }

    //  write_data_to_file("__relxill_unit_returnrad-bbody-xillspec.dat",
    //                    xillspec->ener, xillspec->flu[1], xillspec->n_ener);

    free_xill_spec(xillspec);
    //  free_all_cached_xillTables();
}


TEST_CASE(" Test relxillBB-xill Model ", "[rrad-bbody]")
{
    LocalModel lmod(ModelName::relxillBBxill);
    DefaultSpec default_spec{};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();
    REQUIRE_NOTHROW(lmod.eval_model(spec));
}
