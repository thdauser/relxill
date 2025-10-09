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

extern "C" {
#include "relutility.h"
#include "writeOutfiles.h"
#include "xilltable.h"
}


#define PREC 1e-6


// ------- //
TEST_CASE(" Test Loading Xillver RTR Table ", "[rrad-bbody]")
{
    int status = EXIT_SUCCESS;

    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.get_xill_params();


    xillTable* xill_tab = nullptr;
    REQUIRE_NOTHROW(get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, &status));
}


