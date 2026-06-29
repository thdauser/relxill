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

#include "writeOutfiles.h"

extern "C" {
#include "relutility.h"
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

TEST_CASE (
" rebin_flux_mean_shifted: flat input stays flat for shift=1 "
,
"[rrad-bbody]"
)
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

TEST_CASE (
" rebin_flux_mean_shifted: identity mapping for shift=1 "
,
"[rrad-bbody]"
)
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

TEST_CASE (
" rebin_flux_mean_shifted: large shift zeroes all output bins "
,
"[rrad-bbody]"
)
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

TEST_CASE (
" rebin_flux_mean_shifted: moderate shift zeros high-energy output bins "
,
"[rrad-bbody]"
)
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


TEST_CASE (
" rebin_flux_mean_shifted: matches rebin_mean_flux with explicit shifted grid "
,
"[rrad-bbody]"
)
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

static double integ_flux(const double *ener, int nener, const double *flux) {
  double sum = 0.0;
  for (int ii = 0; ii < nener; ii++) {
    sum += flux[ii] * (ener[ii + 1] - ener[ii]);
  }
  return sum;
}

TEST_CASE (
" check of manually shfting the bbody spectrum gives the same result as if shifting the energy grid "
,
"[rrad-bbody]"
)
 {
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

TEST_CASE (
" Check that 1st and 2nd order returning radiation (mirrored) spectra are different "
,
"[rrad-bbody-norder]"
)
 {

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


TEST_CASE (
" Check if changing the Temperature Profile give different results"
,
"[rrad-bbody-tprofile]"
)
 {

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


TEST_CASE (
" Test Loading Single Xillver from RTR Table should fail "
,
"[rrad-bbody]"
)
{
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.xill_params();

    // by default, xill_param->radzone == -1 should fail to load the xillver spectrum
    REQUIRE(get_xillver_spectra(xill_param, &status) == nullptr);
}

TEST_CASE(
    " Verify that the number of radial zones from the xillver RTR table matches the radial grid in relxill ",
    "[rrad-bbody]"
) {
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.xill_params();
    auto rel_param = lmod.rel_params();

    xillTable *tab = nullptr;
    const char *fname = get_init_xillver_table(&tab, xill_param->model_type, xill_param->prim_type, &status);
    REQUIRE(status == EXIT_SUCCESS);

    // get a standard grid for the convolution (is rebinned later to the input grid)
    auto ener_grid = get_relxill_conv_energy_grid();
    int n_ener = ener_grid->nbins;
    double *ener = ener_grid->ener;
    returnSpec2D *returnSpec = spec_returnrad_blackbody(ener,
                                                        n_ener,
                                                        xill_param->kTbb,
                                                        rel_param->rin,
                                                        rel_param->rout,
                                                        rel_param->a,
                                                        &status);
    int nrad = returnSpec->nrad;


    // test the number of radial zones the model wants is identical to the xillver table parameter values
    for (int ii = 0; ii < tab->num_param; ii++) {
        if (tab->param_index[ii] == PARAM_RID) {
            REQUIRE(tab->num_param_vals[ii] == nrad);
        }
    }
}


// ------- //
TEST_CASE (
" Test Loading Xillver RTR Table "
,
"[rrad-bbody]"
)
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


TEST_CASE (
" Test relxillBB-xill Model "
,
"[relxillbbxill]"
)
{
    LocalModel lmod(ModelName::relxillBBxill);
    DefaultSpec default_spec{};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();
    REQUIRE_NOTHROW(lmod.eval_model(spec));
}

TEST_CASE(
    " Test xillver_RTR spectra normalization ",
    "[xillverRTR]"
) {
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBBxill);
    auto xill_param = lmod.xill_params();

    const int nzone = 50;
    xillParam xpar = *xill_param;
    for (int ii = 0; ii < nzone; ii++) {
        xpar.radzone = ii;
        auto spec = get_xillver_spectra(&xpar, &status);
        double norm_fac = get_avg_xillver_spec_normalization(spec, &status);

        printf(" %i: %.3e \n", ii, norm_fac);
    }
}


// calc_angle_averaged_xillver_flux() has external linkage but is not exposed in
// Xillspec.h; forward-declare it here so the equivalence test can call it.
double *calc_angle_averaged_xillver_flux(const xillSpec *xill_spec, int *status);

TEST_CASE(
    " angle-averaged xillver flux equals angle-weighted intensity with slab flux weights ",
    "[xillverRTR][angledist]"
) {
    // Purpose: prove that calc_angle_averaged_xillver_flux() is exactly the special
    // case of the angle-weighting machinery (calc_xillver_angdep /
    // calculateAngleWeightedXillverSpectrum) when the weights are set to the
    // semi-infinite-slab flux weights  w_i = 0.5*cos(incl_i) / n_incl, which carry
    // the intrinsic-intensity -> observed-flux conversion.
    //
    // It also pins down the (intentional) fact that these slab flux weights do NOT
    // sum to 1 -- unlike rel_cosne->dist, which is a normalized probability distribution.
    int status = EXIT_SUCCESS;

    // a relxill-type model => get_xillver_spectra returns spectra for ALL table inclinations
    // use large inner radius to check for the convergence at large radii, in this case the transfer function
    // should converge towards F=0.5*cos(incl) (see Speith Diploma Thesis, Relxill Internal Docu)
    LocalModel lmod(ModelName::relxill);
    lmod.set_par(XPar::rin, 200);
    lmod.set_par(XPar::rout, 1000);
    xillSpec *spec = get_xillver_spectra(lmod.xill_params(), &status);
    REQUIRE(status == EXIT_SUCCESS);
    REQUIRE(spec != nullptr);
    REQUIRE(spec->n_incl > 1); // averaging is non-trivial only with several inclinations

    // (1) reference: the dedicated angle-averaged-flux routine
    double *spec_avg = calc_angle_averaged_xillver_flux(spec, &status);
    REQUIRE(status == EXIT_SUCCESS);
    REQUIRE(spec_avg != nullptr);

    // build the slab flux weights that reproduce calc_angle_averaged_xillver_flux:
    //   w_i = norm_factor_semi_infinite_slab(incl_i[deg]) / n_incl
    VecD dist(spec->n_incl);
    double dist_sum = 0.0;
    for (int ii = 0; ii < spec->n_incl; ii++) {
        dist[ii] = norm_factor_semi_infinite_slab(spec->incl[ii] * 180.0 / M_PI) / spec->n_incl;
        dist_sum += dist[ii];
    }

    // these slab weights are deliberately NOT a normalized probability distribution
    REQUIRE(fabs(dist_sum - 1.0) > 0.1);

    // (2a) core summation: calc_xillver_angdep with the slab flux weights must reproduce
    //      calc_angle_averaged_xillver_flux bin-by-bin (identical operations => exact)
    VecD spec_weighted(spec->n_ener, 0.0);
    calc_xillver_angdep(spec_weighted.data(), spec, dist.data(), &status);
    REQUIRE(status == EXIT_SUCCESS);

    for (int jj = 0; jj < spec->n_ener; jj++) {
        INFO("calc_xillver_angdep mismatch at bin " << jj);
        REQUIRE(fabs(spec_weighted[jj] - spec_avg[jj]) <= 1e-12 * fabs(spec_avg[jj]));
    }

    // (2b) full public path: calculateAngleWeightedXillverSpectrum rebinned back onto
    //      the native grid (identity rebin) must also reproduce the averaged spectrum
    VecD spec_weighted_rebinned(spec->n_ener, 0.0);
    calc_angle_weighted_xillver_spec(spec_weighted_rebinned,
                                     spec->ener,
                                     spec->n_ener,
                                     dist.data(),
                                     spec,
                                     &status);
    REQUIRE(status == EXIT_SUCCESS);

    for (int jj = 0; jj < spec->n_ener; jj++) {
        INFO("calculateAngleWeightedXillverSpectrum mismatch at bin " << jj);
        REQUIRE(fabs(spec_weighted_rebinned[jj] - spec_avg[jj]) <= 1e-9 * fabs(spec_avg[jj]));
    }

    delete[] spec_avg;
    free_xill_spec(spec);
}


TEST_CASE (
" Test relxillBB Model "
,
"[relxillbb]"
)
{
    LocalModel lmod(ModelName::relxillBB);
    setenv("RELXILL_BBRET_NOREFL", "1", 1);
    DefaultSpec default_spec{};
    XspecSpectrum spec = default_spec.get_xspec_spectrum();
    REQUIRE_NOTHROW(lmod.eval_model(spec));
    unsetenv("RELXILL_BBRET_NOREFL");
}

TEST_CASE (
" Emissivity profile from return radiation spectra "
,
"[rrad-bbody-emisprofile]"
)
 {
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

    // profile should be in ascending radial order
    REQUIRE(emis_profile->re[0] < emis_profile->re[1]);

    // emissivity should be
    // (1) positive and non-zero for all radial zones
    // (2) decreasing as due to returning radiation
    for (int ii = 0; ii < emis_profile->nr; ii++) {
        REQUIRE(emis_profile->emis[ii] > 0.0);
        if (ii > 0) {
            REQUIRE(emis_profile->emis[ii] <= emis_profile->emis[ii-1]);
        }
    }

    free_emisProfile(emis_profile);
    delete_returnSpec2D(&returnSpec);
 }


TEST_CASE(
    " Ionization gradient from return radiation has max equal to input logxi ",
    "[rrad-bbody-iongradient]"
) {
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

    const double input_logxi = 3.0;
    LocalModel lmod(ModelName::relxillBB);
    lmod.set_par(XPar::logxi, input_logxi);

    auto ion_gradient = calc_ion_gradient_from_returnrad(returnSpec, lmod.get_model_params(), &status);
    REQUIRE(status == EXIT_SUCCESS);

    double lxi_max = ion_gradient.lxi[0];
    for (int ii = 1; ii < ion_gradient.nzones(); ii++) {
        if (ion_gradient.lxi[ii] > lxi_max) {
            lxi_max = ion_gradient.lxi[ii];
        }
    }

    REQUIRE(lxi_max == Catch::Approx(input_logxi).epsilon(1e-2));

    delete_returnSpec2D(&returnSpec);
}