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
#include "Rebin.h"
#include "Relreturn_BlackBody.h"
#include "Rellp.h"
#include "Relxill_errorhandling.h"

extern "C" {
#include "relutility.h"
#include "writeOutfiles.h"
#include "xilltable_reading.h"
}


#define PREC 1e-6


// Helper: uniform energy grid with nener bins over [emin, emax]
static VecD make_ener_grid(int nener, double emin, double emax) {
    VecD ener(nener + 1);
    const double de = (emax - emin) / nener;
    for (int i = 0; i <= nener; i++) { ener[i] = emin + i * de; }
    return ener;
}

TEST_CASE(" rebin_flux_mean_shifted: flat input stays flat for shift=1 ", "[rrad-bbody]")
{
    const int nener = 50;
    VecD ener = make_ener_grid(nener, 1.0, 6.0);

    const double C = 3.14;
    VecD flux0(nener, C);
    VecD flux_out(nener, 0.0);

    rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 1.0);

    // bin 0 is 0: known edge effect of the midpoint interpolation scheme
    // (output midpoint coincides exactly with the lower boundary of the input midpoint bracket)
    REQUIRE(fabs(flux_out[0]) < PREC);
    for (int in = 1; in < nener; in++) {
        INFO("Failing index: " << in);
        REQUIRE(fabs(flux_out[in] - C) < PREC);
    }
}

TEST_CASE(" rebin_flux_mean_shifted: identity mapping for shift=1 ", "[rrad-bbody]")
{
    const int nener = 50;
    VecD ener = make_ener_grid(nener, 1.0, 6.0);

    VecD flux0(nener);
    for (int i = 0; i < nener; i++) { flux0[i] = (i + 1) * 0.1; }

    VecD flux_out(nener, 0.0);
    rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 1.0);

    for (int in = 1; in < nener; in++) {
        INFO("Failing index: " << in);
        REQUIRE(fabs(flux_out[in] - flux0[in]) < PREC);
    }
}

TEST_CASE(" rebin_flux_mean_shifted: large shift zeroes all output bins ", "[rrad-bbody]")
{
    const int nener = 50;
    VecD ener = make_ener_grid(nener, 1.0, 6.0);

    VecD flux0(nener, 1.0);
    VecD flux_out(nener, 0.0);

    // shift >> 1: the output grid ener*shift lies entirely above the input grid
    // -> exercises the early-return path; all output bins must be zero
    rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 100.0);

    for (int in = 0; in < nener; in++) {
        INFO("Failing index: " << in);
        REQUIRE(fabs(flux_out[in]) < PREC);
    }
}

TEST_CASE(" rebin_flux_mean_shifted: moderate shift zeros high-energy output bins ", "[rrad-bbody]")
{
    const int nener = 50;
    VecD ener = make_ener_grid(nener, 1.0, 6.0);

    VecD flux0(nener, 1.0);
    VecD flux_out(nener, 0.0);

    // shift=1.5: low-index output bins fall within the input grid (non-zero),
    // high-index bins exceed it (zero via early return)
    rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, 1.5);

    REQUIRE(fabs(flux_out[0] - 1.0) < PREC);         // in range, flat flux preserved
    REQUIRE(fabs(flux_out[nener - 1]) < PREC);        // out of range
}


TEST_CASE(" rebin_flux_mean_shifted: matches rebin_mean_flux with explicit shifted grid ", "[rrad-bbody]")
{
    const int nener = 50;
    VecD ener = make_ener_grid(nener, 1.0, 6.0);

    // arbitrary non-trivial flux values
    VecD flux0(nener);
    for (int i = 0; i < nener; i++) { flux0[i] = 1.0 + 0.1 * i + 0.002 * i * i; }

    const double shift = 0.8;

    // manually build the shifted output grid
    VecD ener_shifted(nener + 1);
    for (int i = 0; i <= nener; i++) { ener_shifted[i] = ener[i] * shift; }

    // reference: rebin_mean_flux with the explicitly shifted grid as output
    VecD flux_ref(nener, 0.0);
    int status = EXIT_SUCCESS;
    rebin_mean_flux(ener_shifted.data(), flux_ref.data(), nener,
                    ener.data(), flux0.data(), nener, &status);
    REQUIRE(status == EXIT_SUCCESS);

    // function under test
    VecD flux_out(nener, 0.0);
    rebin_mean_flux_shifted(ener.data(), nener, flux_out, flux0, shift);

    for (int in = 0; in < nener; in++) {
        INFO("Failing index: " << in);
        REQUIRE(fabs(flux_out[in] - flux_ref[in]) < PREC);
    }
}

static double integ_flux(const double* ener, int nener, const double* flux) {
    double sum = 0.0;
    for (int ii = 0; ii < nener; ii++) {
        sum += flux[ii] * (ener[ii + 1] - ener[ii]);
    }
    return sum;
}

TEST_CASE(" check of manually shfting the bbody spectrum gives the same result as if shifting the energy grid ", "[rrad-bbody]") {
    const int nener = 5000;

    double ener[nener + 1];
    get_log_grid(ener, nener, 0.01, 150.0);

    double gfac = 1.5;
    double temperature = 1.0;
    VecD spec(nener);
    bbody_spec(ener,  nener, spec.data(), temperature, gfac);

    VecD spec0(nener);
    bbody_spec(ener,  nener, spec0.data(), temperature, 1.0);
    VecD spec0_shift(nener);
    // as we shift the energy grid and not the spectrum, we need to shift it backwards with 1/gfac
    rebin_mean_flux_shifted(ener, nener, spec0_shift, spec0, 1./gfac);

    auto pho_flux = integ_flux(ener, nener, spec.data());
    auto pho_flux0 = integ_flux(ener, nener, spec0_shift.data());

    REQUIRE(fabs(pho_flux-pho_flux0) < 1e-3);

}

TEST_CASE(" Check that 1st and 2nd order returning radiation (mirrored) spectra are different ", "[rrad-bbody-norder]") {

    LocalModel lmod(ModelName::relxillBBxill);
    DefaultSpec default_spec{};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();

    setenv("RELXILL_BBRET_NOREFL","1",1);
    setenv("RELXILL_RRAD_BBODY_ORDER","1",1);

    lmod.eval_model(spec);
    double sum_1st = calcSum(spec.flux, spec.num_flux_bins());

    setenv("RELXILL_RRAD_BBODY_ORDER","2",1);
    lmod.eval_model(spec);
    double sum_2nd = calcSum(spec.flux, spec.num_flux_bins());

    REQUIRE((sum_2nd - sum_1st) > 1e-3);

    unsetenv("RELXILL_BBRET_NOREFL");
    unsetenv("RELXILL_RRAD_BBODY_ORDER");
}


TEST_CASE(" Check if changing the Temperature Profile give different results", "[rrad-bbody-tprofile]") {

    LocalModel lmod(ModelName::relxillBBxill);
    DefaultSpec default_spec{};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();

    setenv("RELXILL_BBRET_NOREFL","1",1);

    lmod.eval_model(spec);
    double sum_default = calcSum(spec.flux, spec.num_flux_bins());

    setenv("RELXILL_RRAD_BBODY_TPROFILE","1",1);
    lmod.eval_model(spec);
    double sum_1 = calcSum(spec.flux, spec.num_flux_bins());

    // check of setting it to 1 is the default
    REQUIRE(fabs(sum_default - sum_1) < 1e-6);


    setenv("RELXILL_RRAD_BBODY_TPROFILE","2",1);
    lmod.eval_model(spec);
    double sum_2 = calcSum(spec.flux, spec.num_flux_bins());

    unsetenv("RELXILL_BBRET_NOREFL");
    unsetenv("RELXILL_RRAD_BBODY_TPROFILE");

    // check that 1 and 2 are different
    REQUIRE(fabs(sum_2 - sum_1) > 1e-6);
}


TEST_CASE(" Test Loading Single Xillver from RTR Table should fail ", "[rrad-bbody]")
{
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.xill_params();

    // by default, xill_param->radzone == -1 should fail to load the xillver spectrum
    REQUIRE(get_xillver_spectra(xill_param, &status) == nullptr);
}


// ------- //
TEST_CASE(" Test Loading Xillver RTR Table ", "[rrad-bbody]")
{
    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.get_xill_param_struct();

    int status = EXIT_SUCCESS;
    xillTable* xill_tab = nullptr;
    REQUIRE_NOTHROW([&]()
        {
        auto tmp = get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, &status);
        }()
    );

    xillSpec* xillspec = nullptr;
    xill_param->radzone = 2;
    xillspec = get_xillver_spectra(xill_param, &status);
    REQUIRE(xillspec != nullptr);

     for (int ii = 0; ii < xillspec->n_incl; ii++)
    {
        INFO("Failing index: " << ii);
        REQUIRE(sum_flux(xillspec->flu[ii], xillspec->n_ener) > 1e-10);
    }

      write_data_to_file("__relxill_unit_returnrad-bbody-xillspec.dat",
                         xillspec->ener, xillspec->flu[1], xillspec->n_ener);

    delete_xill_params(xill_param);
    free_xill_spec(xillspec);
    //  free_all_cached_xillTables();
}


TEST_CASE(" Test relxillBB-xill Model ", "[relxillbbxill]")
{
    LocalModel lmod(ModelName::relxillBBxill);
    DefaultSpec default_spec{};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();
    REQUIRE_NOTHROW(lmod.eval_model(spec));
}

TEST_CASE(" Test relxillBB Model ", "[relxillbb]")
{
    LocalModel lmod(ModelName::relxillBB);
    setenv("RELXILL_BBRET_NOREFL", "1", 1);
    DefaultSpec default_spec{};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();
    REQUIRE_NOTHROW(lmod.eval_model(spec));
    unsetenv("RELXILL_BBRET_NOREFL");
}

TEST_CASE(" Emissivity profile from return radiation spectra ", "[rrad-bbody-emisprofile]") {
    int status = EXIT_SUCCESS;

    double spin = 0.998;
    double Rin = kerr_rms(spin);
    double Rout = RMAX_RELRET;
    double Tin = 2.0;

    int n_ener = 500;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.01, 50.0);

    returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, n_ener, Tin, Rin, Rout, spin, &status);
    REQUIRE(status == EXIT_SUCCESS);
    REQUIRE(returnSpec != nullptr);

    emisProfile *emis_profile =
            calculate_emis_profile_from_spectra(returnSpec->specRet,
                                                returnSpec->rlo,
                                                returnSpec->rhi,
                                                returnSpec->nrad,
                                                returnSpec->ener,
                                                returnSpec->n_ener,
                                                &status);
    REQUIRE(status == EXIT_SUCCESS);
    REQUIRE(emis_profile != nullptr);
    REQUIRE(emis_profile->nr == returnSpec->nrad);

    // emissivity should be positive and non-zero for all radial zones
    for (int ii = 0; ii < emis_profile->nr; ii++) {
        REQUIRE(emis_profile->emis[ii] > 0.0);
    }

    free_emisProfile(emis_profile);
    delete_returnSpec2D(&returnSpec);
}
