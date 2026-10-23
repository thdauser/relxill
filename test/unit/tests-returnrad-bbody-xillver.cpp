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


// Unit tests of the reflection of the returning thermal radiation (relxillBB, relxillBB_NS): loading of the
// xillverRTR (radial zones) and xillverNS tables, and the angle weighting and normalization of the xillver spectra of
// the zones (see doc/returnrad_normalizations.md). Run with bin/tests "[rrad-bbody][xillver]".

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "Xillspec.h"
#include "Xilltable_handling.h"
#include "Relreturn_BlackBody.h"
#include "Relbase.h"
#include "Relprofile.h"
#include "Rebin.h"
#include "common-functions.h"
#include "writeOutfiles.h"

extern "C" {
#include "relutility.h"
#include "xilltable_reading.h"
}

// xillver normalization of a spectrum, computed as for N_xill in
// calc_zone_xillver_reflected_flux(): rebin onto the coarse xillver grid,
// then calc_xillver_norm_from_primary_spectrum().
static double get_xillver_norm_on_coarse_grid(const double *spec, const double *ener, int n_ener, int *status) {
    EnerGrid *egrid = get_coarse_xillver_energrid();
    VecD spec_coarse(egrid->nbins, 0.0);
    _rebin_spectrum(egrid->ener, spec_coarse.data(), egrid->nbins, ener, spec, n_ener);
    return calc_xillver_norm_from_primary_spectrum(spec_coarse.data(), egrid->ener, egrid->nbins, status);
}

// per-zone ratio N(reflect)/N(return) of the output of calc_zone_xillver_reflected_flux().
// boost=-1: |boost|=1 does not scale the reflection
//  => xill_flux_returnrad = N_xill(specRet) * F_xill, with the angle weighting given by rel_profile->rel_cosne
static VecD calc_ratio_reflected_to_returning_flux(const xillParam &xill_param,
                                                   const relline_spec_multizone *rel_profile,
                                                   const returnSpec2D *returnSpec,
                                                   int *status) {
    xillParam xill_param_zone = xill_param;
    xill_param_zone.boost = -1.0;

    VecD ratio(returnSpec->nrad, 0.0);
    VecD xill_flux_returnrad(returnSpec->n_ener, 0.0);

    // printf("\n %5s  %12s  %12s  %12s\n", "izone", "N(return)", "N(reflect)", "refl/return");
    for (int izone = 0; izone < returnSpec->nrad; izone++) {
        xill_param_zone.radzone = static_cast<double>(izone);
        calc_zone_xillver_reflected_flux(xill_param_zone,
                                         rel_profile,
                                         returnSpec,
                                         xill_flux_returnrad,
                                         izone,
                                         status);

        double norm_return =
                get_xillver_norm_on_coarse_grid(returnSpec->specRet[izone],
                                                returnSpec->ener,
                                                returnSpec->n_ener,
                                                status);
        double norm_reflect =
                get_xillver_norm_on_coarse_grid(xill_flux_returnrad.data(),
                                                returnSpec->ener,
                                                returnSpec->n_ener,
                                                status);
        if (*status != EXIT_SUCCESS) {
            break;
        }

        ratio[izone] = norm_reflect / norm_return;
        //     printf(" %5i  %12.4e  %12.4e  %12.4e\n", izone, norm_return, norm_reflect, ratio[izone]);
    }
    return ratio;
}

// relline_spec_multizone with only the emission-angle distribution rel_cosne set (the only part
// used by calc_zone_xillver_reflected_flux()), owning the memory of dist [n_zones][n_cosne]
struct AngleDistProfile {
    VecD cosne;
    VecD_2D dist;
    std::vector<double *> dist_rows;
    RelCosne rel_cosne{};
    relline_spec_multizone profile{};

    AngleDistProfile(const VecD &cosne_in, int n_zones)
        : cosne(cosne_in), dist(n_zones, VecD(cosne_in.size(), 0.0)), dist_rows(n_zones) {
        for (int ii = 0; ii < n_zones; ii++) {
            dist_rows[ii] = dist[ii].data();
        }
        rel_cosne.n_cosne = static_cast<int>(cosne.size());
        rel_cosne.n_zones = n_zones;
        rel_cosne.cosne = cosne.data();
        rel_cosne.dist = dist_rows.data();
        profile.n_zones = n_zones;
        profile.rel_cosne = &rel_cosne;
    }

    // holds pointers to its own members
    AngleDistProfile(const AngleDistProfile &) = delete;

    AngleDistProfile &operator=(const AngleDistProfile &) = delete;
};

// semi-infinite slab flux weights 0.5*cos(incl_i)/n_incl of the table inclinations (xillSpec::incl is
// given in deg), computed independently of get_semi_infinite_slab_flux_weights() to test it
static VecD expected_slab_flux_weights(const xillSpec *xill_spec) {
    VecD weights(xill_spec->n_incl);
    for (int ii = 0; ii < xill_spec->n_incl; ii++) {
        weights[ii] = 0.5 * cos(xill_spec->incl[ii] * M_PI / 180.0) / xill_spec->n_incl;
    }
    return weights;
}

TEST_CASE (
"xillverRTR table: loading spectra without a radial zone (radzone=-1) fails"
,
"[rrad-bbody][xillver][xillverRTR]"
)
 {
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBB);
    auto xill_param = lmod.xill_params();

    // by default, xill_param->radzone == -1 should fail to load the xillver spectrum
    REQUIRE(get_xillver_spectra(xill_param, &status) == nullptr);
}

TEST_CASE (
"xillverRTR table: number of radial zones equals the radial grid of the returning radiation"
,
"[rrad-bbody][xillver][xillverRTR]"
)
 {
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBB);
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

    delete_returnSpec2D(&returnSpec);
}

TEST_CASE (
"xillverRTR table: spectra of a radial zone are loaded and non-zero"
,
"[rrad-bbody][xillver][xillverRTR]"
)
 {
    LocalModel lmod(ModelName::relxillBB);
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

    for (int ii = 0; ii < xillspec->n_incl; ii++) {
        INFO("Failing index: " << ii);
        REQUIRE(sum_flux(xillspec->flu[ii], xillspec->n_ener) > 1e-10);
    }

    write_data_to_file("__relxill_unit_returnrad-bbody-xillspec.dat",
                       xillspec->ener,
                       xillspec->flu[1],
                       xillspec->n_ener);

    delete_xill_params(xill_param);
    free_xill_spec(xillspec);
    //  free_all_cached_xillTables();
}

TEST_CASE (
"xillverNS table: spectra for relxillBB_NS are loaded and non-zero"
,
"[rrad-bbody][xillver]"
)
 {
  int status = EXIT_SUCCESS;

  LocalModel lmod(ModelName::relxillBB_NS);
  auto xill_params = lmod.xill_params();
  auto xill_table_params = get_xilltab_param(xill_params, &status);
  auto xill_spec = get_xillver_spectra_table(xill_table_params, &status);

  for (size_t ii = 0; ii < xill_spec->n_incl; ii++)
  {
    REQUIRE(calcSum(xill_spec->flu[ii], xill_spec->n_ener) > 1e-8);
  }

  free_xill_spec(xill_spec);
  free(xill_table_params);
}

// ------- //
// Semi-infinite slab flux weights 0.5*cos(incl_i)/n_incl of the xillver table inclinations, which convert the
// intrinsic xillver intensities into the angle-averaged flux: get_semi_infinite_slab_flux_weights() and
// calc_angle_averaged_xillver_flux() (used for the returning radiation correction factors of the lamp post models)
// ------- //
TEST_CASE (
"xillver angle weighting: semi-infinite slab flux weights 0.5*cos(incl)/n_incl and the angle-averaged flux"
,
"[rrad-bbody][xillver][xillverRTR][slabweights]"
)
 {
    // Tests get_semi_infinite_slab_flux_weights() and calc_angle_averaged_xillver_flux() against independently
    // computed values. As the table inclinations are the bin centres of a grid uniform in cos(incl), the weights sum
    // to int_0^1 0.5*mu dmu = 0.25. (The weights were wrong before, as the table inclinations, already given in deg,
    // were converted to deg again.)
    int status = EXIT_SUCCESS;

    // relxill: standard xillver table; relxillBB: xillverRTR table
    const ModelName mod_name = GENERATE(ModelName::relxill, ModelName::relxillBB);
    LocalModel lmod(mod_name);
    INFO("model " << lmod.get_model_string());

    xillParam *xill_param = lmod.get_xill_param_struct();
    xill_param->radzone = 0; // radial zone of the xillverRTR table (not used by the standard table)
    xillSpec *xill_spec = get_xillver_spectra(xill_param, &status);
    REQUIRE(status == EXIT_SUCCESS);
    REQUIRE(xill_spec->n_incl > 1);

    const VecD weights = get_semi_infinite_slab_flux_weights(xill_spec);
    const VecD weights_expected = expected_slab_flux_weights(xill_spec);
    REQUIRE(weights.size() == weights_expected.size());

    double sum_weights = 0.0;
    for (int ii = 0; ii < xill_spec->n_incl; ii++) {
        INFO("table inclination " << xill_spec->incl[ii] << " deg");
        CHECK(weights[ii] == Catch::Approx(weights_expected[ii]).epsilon(1e-10));
        sum_weights += weights[ii];
    }
    CHECK(sum_weights == Catch::Approx(0.25).epsilon(1e-3));

    // the angle-averaged flux is the sum of the table spectra weighted with these weights
    double *spec_avg = calc_angle_averaged_xillver_flux(xill_spec, &status);
    REQUIRE(status == EXIT_SUCCESS);
    for (int jj = 0; jj < xill_spec->n_ener; jj++) {
        double spec_avg_expected = 0.0;
        for (int ii = 0; ii < xill_spec->n_incl; ii++) {
            spec_avg_expected += weights_expected[ii] * xill_spec->flu[ii][jj];
        }
        INFO("energy bin " << jj);
        CHECK(spec_avg[jj] == Catch::Approx(spec_avg_expected).epsilon(1e-10));
    }

    delete[] spec_avg;
    free_xill_spec(xill_spec);
    delete_xill_params(xill_param);
}

// ------- //
// relxillBB angle weighting of the xillver_RTR spectra: normalization of the GR emission-angle distribution,
// energy conservation of the reflected spectrum (angle-averaged flux vs. returning flux), and comparison of the
// angle weightings (angle-mixed intensity vs. angle-averaged flux)
// ------- //
TEST_CASE (
"xillverRTR angle weighting: GR emission-angle distribution of each zone sums to 1"
,
"[rrad-bbody][xillver][xillverRTR][cosnedist]"
)
 {
    // The emission-angle distribution weights the xillver spectra of the table inclinations in each radial zone
    // (calc_xillver_angdep()) and is a probability distribution (dist >= 0, sum = 1). Only the inclination grid of
    // the xillver table and the radial zones of the returning radiation are needed, no xillver spectra.
    int status = EXIT_SUCCESS;
    const double incl_deg = GENERATE(5.0, 60.0);
    INFO("incl = " << incl_deg << " deg");

    LocalModel lmod(ModelName::relxillBB);
    lmod.set_par(XPar::incl, incl_deg);
    relParam *rel_param = lmod.get_rel_params_struct();
    xillParam *xill_param = lmod.get_xill_param_struct();

    // same set-up as in relxillbb_kernel()
    auto ener_grid = get_relxill_conv_energy_grid();
    int n_ener = ener_grid->nbins;
    double *ener = ener_grid->ener;

    xillTable *xill_tab = nullptr;
    get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, &status);
    REQUIRE(status == EXIT_SUCCESS);

    // the radial zones are the ones of the returning radiation
    returnSpec2D *returnSpec = spec_returnrad_blackbody(ener,
                                                        n_ener,
                                                        xill_param->kTbb,
                                                        rel_param->rin,
                                                        rel_param->rout,
                                                        rel_param->a,
                                                        &status);
    REQUIRE(status == EXIT_SUCCESS);
    auto radialGrid = get_radial_grid_from_returntab(returnSpec);

    RelSysPar *systemParameters = get_system_parameters(rel_param, &status);
    relline_spec_multizone *relProfile = relbase_profile(ener,
                                                         n_ener,
                                                         rel_param,
                                                         systemParameters,
                                                         xill_tab,
                                                         radialGrid,
                                                         &status);
    REQUIRE(status == EXIT_SUCCESS);
    REQUIRE(relProfile->n_zones == returnSpec->nrad);

    const RelCosne *rel_cosne = relProfile->rel_cosne;
    REQUIRE(rel_cosne != nullptr);
    REQUIRE(rel_cosne->n_cosne == xill_tab->n_incl);

    // the emission-angle distribution is a probability distribution in each zone
    for (int izone = 0; izone < rel_cosne->n_zones; izone++) {
        double sum = 0.0;
        for (int ii = 0; ii < rel_cosne->n_cosne; ii++) {
            INFO("zone " << izone << ", cosne bin " << ii);
            CHECK(rel_cosne->dist[izone][ii] >= 0.0);
            sum += rel_cosne->dist[izone][ii];
        }
        INFO("zone " << izone);
        CHECK(sum == Catch::Approx(1.0).epsilon(1e-8));
    }

    delete_returnSpec2D(&returnSpec);
    delete_rel_params(rel_param);
    delete_xill_params(xill_param);
}

TEST_CASE (
"xillverRTR normalization: reflected flux is below the returning flux in each zone"
,
"[rrad-bbody][xillver][xillverRTR][norm]"
)
 {
    // Energy conservation: the GR emission-angle distribution is replaced by the semi-infinite slab flux weights
    // (identical for all zones), which turns the angle-weighted intensity into the angle-averaged flux (as in
    // calc_angle_averaged_xillver_flux()). This reflected flux must be below the returning flux. (With the GR
    // emission-angle distribution, the result is the angle-mixed intensity, which can not be compared with the
    // returning flux, see "[angleweight]".)
    int status = EXIT_SUCCESS;
    LocalModel lmod(ModelName::relxillBB);
    relParam *rel_param = lmod.get_rel_params_struct();
    xillParam *xill_param = lmod.get_xill_param_struct();

    // same set-up as in relxillbb_kernel()
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
    REQUIRE(status == EXIT_SUCCESS);

    xillParam xill_param_incl = *xill_param;
    xill_param_incl.radzone = 0;
    xillSpec *xill_spec = get_xillver_spectra(&xill_param_incl, &status);
    REQUIRE(status == EXIT_SUCCESS);

    // slab flux weights 0.5*cos(incl_i)/n_incl (the weights themselves are tested in "[slabweights]")
    const VecD slab_weights = get_semi_infinite_slab_flux_weights(xill_spec);

    VecD cosne(xill_spec->n_incl);
    for (int ii = 0; ii < xill_spec->n_incl; ii++) {
        cosne[ii] = cos(xill_spec->incl[ii] * M_PI / 180.0);
    }
    free_xill_spec(xill_spec);

    AngleDistProfile slab(cosne, returnSpec->nrad);
    for (int izone = 0; izone < returnSpec->nrad; izone++) {
        std::copy(slab_weights.begin(), slab_weights.end(), slab.dist[izone].begin()); // keep dist_rows valid
    }

    VecD ratio = calc_ratio_reflected_to_returning_flux(*xill_param, &slab.profile, returnSpec, &status);
    REQUIRE(status == EXIT_SUCCESS);

    // energy conservation: the reflected flux is below the returning (illuminating) flux
    for (int izone = 0; izone < returnSpec->nrad; izone++) {
        INFO("zone " << izone);
        CHECK(ratio[izone] > 0.0);
        CHECK(ratio[izone] < 1.0);
    }

    delete_returnSpec2D(&returnSpec);
    delete_rel_params(rel_param);
    delete_xill_params(xill_param);
}

TEST_CASE (
"xillverRTR angle weighting: angle-mixed intensity vs. angle-averaged flux per zone (written to file)"
,
"[rrad-bbody][xillver][xillverRTR][angleweight]"
)
 {
    // The reflected spectrum per zone (as in rframe-reflect.fits of the diagnose outfiles) is the angle-mixed
    // intensity, weighted with the GR emission-angle distribution dist of the observer inclination. Only the
    // angle-averaged flux of the semi-infinite slab can be compared with the returning flux (energy conservation).
    // With r_i = N(I_i) the xillver normalization of the table spectrum at inclination i (as reflect =
    // N_xill(ret) * F_xill, it is N(reflect)/N(return) = N(F_xill)), the following weightings are compared:
    //   A_mix  = sum_i dist_i r_i                    angle-mixed intensity (what is convolved in the model)
    //   B_slab = sum_i 0.5*cos(incl_i)/n_incl r_i    angle-averaged flux of the semi-infinite slab
    //   C_obs  = 0.5*cos(incl_obs) * A_mix           slab factor of the observer inclination
    //   D_GR   = sum_i dist_i 0.5*cos(incl_i) r_i    slab factor per emission angle
    // The values per zone are written to the file below (sums are weighted with the returning flux N_ret).
    const char *fname = "__relxill_unit_xillverRTR-angleweight.dat";
    FILE *fp = fopen(fname, "w");
    REQUIRE(fp != nullptr);
    fprintf(fp, "# incl[deg] zone r_mid[Rg] N_ret mu_GR A_mix B_slab C_obs D_GR r_i(table incl, ascending)\n");

    int status = EXIT_SUCCESS;
    for (double incl_deg: {20.0, 40.0, 60.0, 80.0}) {
        INFO("incl = " << incl_deg << " deg");
        LocalModel lmod(ModelName::relxillBB);
        lmod.set_par(XPar::incl, incl_deg);
        relParam *rel_param = lmod.get_rel_params_struct();
        xillParam *xill_param = lmod.get_xill_param_struct();

        // same set-up as in relxillbb_kernel()
        auto ener_grid = get_relxill_conv_energy_grid();
        int n_ener = ener_grid->nbins;
        double *ener = ener_grid->ener;

        xillTable *xill_tab = nullptr;
        get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, &status);
        returnSpec2D *returnSpec = spec_returnrad_blackbody(ener,
                                                            n_ener,
                                                            xill_param->kTbb,
                                                            rel_param->rin,
                                                            rel_param->rout,
                                                            rel_param->a,
                                                            &status);
        REQUIRE(status == EXIT_SUCCESS);

        auto radialGrid = get_radial_grid_from_returntab(returnSpec);
        RelSysPar *systemParameters = get_system_parameters(rel_param, &status);
        relline_spec_multizone *relProfile = relbase_profile(ener,
                                                             n_ener,
                                                             rel_param,
                                                             systemParameters,
                                                             xill_tab,
                                                             radialGrid,
                                                             &status);
        REQUIRE(status == EXIT_SUCCESS);
        const RelCosne *rel_cosne = relProfile->rel_cosne;
        REQUIRE(rel_cosne != nullptr);
        REQUIRE(rel_cosne->n_cosne == xill_tab->n_incl);

        double sum_ret = 0.0, sum_mix = 0.0, sum_slab = 0.0, sum_obs = 0.0, sum_gr = 0.0;
        for (int izone = 0; izone < returnSpec->nrad; izone++) {
            xillParam xill_param_zone = *xill_param;
            xill_param_zone.radzone = static_cast<double>(izone);
            xillSpec *xill_spec = get_xillver_spectra(&xill_param_zone, &status);
            REQUIRE(status == EXIT_SUCCESS);

            const VecD slab_weights = expected_slab_flux_weights(xill_spec);
            VecD r(xill_spec->n_incl);
            double A_mix = 0.0, B_slab = 0.0, D_GR = 0.0, mu_GR = 0.0;
            for (int ii = 0; ii < xill_spec->n_incl; ii++) {
                r[ii] = get_xillver_norm_on_coarse_grid(xill_spec->flu[ii],
                                                        xill_spec->ener,
                                                        xill_spec->n_ener,
                                                        &status);
                const double mu = cos(xill_spec->incl[ii] * M_PI / 180.0);
                const double dist = rel_cosne->dist[izone][ii];
                A_mix += dist * r[ii];
                B_slab += slab_weights[ii] * r[ii];
                D_GR += dist * 0.5 * mu * r[ii];
                mu_GR += dist * mu;
            }
            const double C_obs = 0.5 * cos(incl_deg * M_PI / 180.0) * A_mix;
            const double N_ret = get_xillver_norm_on_coarse_grid(returnSpec->specRet[izone],
                                                                 returnSpec->ener,
                                                                 returnSpec->n_ener,
                                                                 &status);
            REQUIRE(status == EXIT_SUCCESS);

            fprintf(fp,
                    "%5.1f %3i %10.4e %10.4e %6.3f %7.4f %7.4f %7.4f %7.4f",
                    incl_deg,
                    izone,
                    0.5 * (returnSpec->rlo[izone] + returnSpec->rhi[izone]),
                    N_ret,
                    mu_GR,
                    A_mix,
                    B_slab,
                    C_obs,
                    D_GR);
            for (int ii = 0; ii < xill_spec->n_incl; ii++) {
                fprintf(fp, " %6.3f", r[ii]);
            }
            fprintf(fp, "\n");
            free_xill_spec(xill_spec);

            sum_ret += N_ret;
            sum_mix += N_ret * A_mix;
            sum_slab += N_ret * B_slab;
            sum_obs += N_ret * C_obs;
            sum_gr += N_ret * D_GR;

            // energy conservation: the angle-averaged reflected flux is below the returning flux
            INFO("zone " << izone << ": A_mix = " << A_mix << ", B_slab = " << B_slab
                << ", C_obs = " << C_obs << ", D_GR = " << D_GR);
            CHECK(B_slab > 0.0);
            CHECK(B_slab < 1.0);
        }
        fprintf(fp,
                "# incl = %.1f deg, all zones: A_mix = %.4f, B_slab = %.4f, C_obs = %.4f, D_GR = %.4f\n",
                incl_deg,
                sum_mix / sum_ret,
                sum_slab / sum_ret,
                sum_obs / sum_ret,
                sum_gr / sum_ret);

        delete_returnSpec2D(&returnSpec);
        delete_rel_params(rel_param);
        delete_xill_params(xill_param);
    }
    fclose(fp);
}
