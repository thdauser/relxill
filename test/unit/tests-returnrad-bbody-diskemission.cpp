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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/

// Unit tests of the direct disk emission of relxillBB / relxillBB_NS (temperature profiles, radial zones, kernel of
// the direct emission, and the consistency of direct and reflected radiation), partly compared to the Xspec models
// diskbb, ezdiskbb and kerrbb (reference files in test/unit/refdata/). See doc/returnrad_disk_emission.md.
// Run with bin/tests "[diskemission]".

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "XspecSpectrum.h"
#include "Relreturn_BlackBody.h"
#include "Xillspec.h"
#include "Rebin.h"
#include "Relphysics.h"

#include <array>

extern "C" {
#include "relutility.h"
#include "xilltable_reading.h"
#include "envoptions.h"
}

#if !defined(RELXILL_SOURCE_DIR)
#warning: RELXILL_SOURCE_DIR is not defined
#endif

#define REFDATA_DIR "test/unit/refdata/"  // relativ to RELXILL_SOURCE_DIR (set when compiling)


// reference file with the diskbb model of Xspec (column "diskbb") for the given Tin, on the energy grid of the
// disk emission tests below (created with test/unit/refdata/create_refvalues_diskbb.sl)
static std::string get_fname_refdata_diskbb(double Tin) {
  char fname[100];
  snprintf(fname, sizeof(fname), "refvalues_diskbb_Tin%.1f.fits", Tin);
  return fname;
}

// reference file with the ezdiskbb model of Xspec (column "ezdiskbb") for the given Tmax, on the energy grid of the
// disk emission tests below (created with test/unit/refdata/create_refvalues_diskbb.sl)
static std::string get_fname_refdata_ezdiskbb(double Tmax) {
  char fname[100];
  snprintf(fname, sizeof(fname), "refvalues_ezdiskbb_Tmax%.1f.fits", Tmax);
  return fname;
}

// reference file with the kerrbb model of Xspec (column "kerrbb") for the given spin and inclination, with a maximal
// temperature of the disk of 1 keV, on the energy grid of the disk emission tests below (created with
// test/unit/refdata/create_refvalues_kerrbb.sl)
static std::string get_fname_refdata_kerrbb(double spin, double incl) {
  char fname[100];
  snprintf(fname, sizeof(fname), "refvalues_kerrbb_a%.6f_i%.0f.fits", spin, incl);
  return fname;
}

// reads a model spectrum of Xspec (column colname) from a reference file and checks that its energy grid
// (column "elo") is identical to ener
static VecD fits_read_refdata_xspec(const std::string &fname_ref,
                                    const char *colname,
                                    const double *ener,
                                    int n_ener,
                                    int *status) {
  VecD spec_ref(n_ener, 0.0);
  VecD elo(n_ener, 0.0);
  const std::string fullname = std::string(RELXILL_SOURCE_DIR) + "/" + REFDATA_DIR + fname_ref;

  fitsfile *fptr = nullptr;
  if (fits_open_table(&fptr, fullname.c_str(), READONLY, status)) {
    relxill_check_fits_error(status);
    return spec_ref;
  }

  long nrows = 0;
  int colnum_elo = 0;
  int colnum_spec = 0;
  char colname_elo[] = "elo";
  std::string colname_spec(colname);
  fits_get_num_rows(fptr, &nrows, status);
  fits_get_colnum(fptr, CASEINSEN, colname_elo, &colnum_elo, status);
  fits_get_colnum(fptr, CASEINSEN, colname_spec.data(), &colnum_spec, status);
  if (*status == EXIT_SUCCESS && nrows == n_ener) {
    int anynul = 0;
    double nullval = 0.0;
    fits_read_col(fptr, TDOUBLE, colnum_elo, 1, 1, nrows, &nullval, elo.data(), &anynul, status);
    fits_read_col(fptr, TDOUBLE, colnum_spec, 1, 1, nrows, &nullval, spec_ref.data(), &anynul, status);
  }
  relxill_check_fits_error(status);
  int status_close = EXIT_SUCCESS;
  fits_close_file(fptr, &status_close);
  CHECK_STATUS_RET(*status, spec_ref);

  if (nrows != n_ener) {
    printf(" *** error : reference file %s has %li bins, but expected %i bins\n", fname_ref.c_str(), nrows, n_ener);
    *status = EXIT_FAILURE;
    return spec_ref;
  }
  for (int ii = 0; ii < n_ener; ii++) {
    if (fabs(elo[ii] / ener[ii] - 1) > 1e-6) {
      printf(" *** error : energy grid of reference file %s differs in bin %i (%e instead of %e keV)\n",
             fname_ref.c_str(),
             ii,
             elo[ii],
             ener[ii]);
      *status = EXIT_FAILURE;
      return spec_ref;
    }
  }
  return spec_ref;
}

// max. relative deviation of the shape of spec from spec_ref in the band [emin, emax], i.e., of the ratio
// spec/spec_ref normalized to its mean in the band (the mean ratio is returned in mean_ratio, if given)
static double max_shape_deviation(const double *ener,
                                  const double *spec,
                                  const VecD &spec_ref,
                                  int n_ener,
                                  double emin,
                                  double emax,
                                  double *mean_ratio = nullptr) {
  VecD ratio;
  for (int ii = 0; ii < n_ener; ii++) {
    if (ener[ii] >= emin && ener[ii + 1] <= emax) {
      ratio.push_back(spec[ii] / spec_ref[ii]);
    }
  }
  REQUIRE(!ratio.empty());

  double mean = 0.0;
  for (double val: ratio) {
    mean += val / static_cast<double>(ratio.size());
  }
  double max_dev = 0.0;
  for (double val: ratio) {
    max_dev = fmax(max_dev, fabs(val / mean - 1));
  }
  if (mean_ratio != nullptr) {
    *mean_ratio = mean;
  }
  return max_dev;
}

// logarithmic radial grid with n zones from rin to rout
static void get_log_radial_grid(double rin, double rout, int n, VecD &rlo, VecD &rhi) {
  rlo.resize(n);
  rhi.resize(n);
  for (int ii = 0; ii < n; ii++) {
    rlo[ii] = rin * pow(rout / rin, static_cast<double>(ii) / n);
    rhi[ii] = rin * pow(rout / rin, static_cast<double>(ii + 1) / n);
  }
}

// multicolor disk spectrum (spec_multicolor_disk) for the radial zones [rlo, rhi] on the grid ener: as in the model,
// it is calculated on the fine energy grid of the convolution and then rebinned (evaluating the blackbody at the bin
// centres of a coarse grid directly would add a deviation of ~3%)
static VecD calc_disk_spectrum_rebinned(const VecD &rlo,
                                        const VecD &rhi,
                                        double kT,
                                        int tprofile,
                                        const double *ener,
                                        int n_ener,
                                        int *status) {
  const EnerGrid *conv_grid = get_relxill_conv_energy_grid();
  VecD spec_conv(conv_grid->nbins);
  spec_multicolor_disk(conv_grid->ener, spec_conv.data(), conv_grid->nbins, rlo, rhi, kT, tprofile, status);
  VecD spec(n_ener, 0.0);
  _rebin_spectrum(ener, spec.data(), n_ener, conv_grid->ener, spec_conv.data(), conv_grid->nbins);
  return spec;
}

// max. shape deviation (see max_shape_deviation) of the direct disk emission (boost=0) of the full relxillBB_NS model
// from spec_ref in the band [emin, emax], for the temperature profile tprofile (value of RELXILL_RRAD_BBODY_TPROFILE)
static double shape_deviation_relxillbb_ns_disk_emission(const char *tprofile,
                                                         double kTbb,
                                                         double spin,
                                                         double rin,
                                                         double incl,
                                                         const double *ener,
                                                         int n_ener,
                                                         const VecD &spec_ref,
                                                         double emin,
                                                         double emax) {
  setenv("RELXILL_RRAD_BBODY_TPROFILE", tprofile, 1);
  LocalModel lmod(ModelName::relxillBB_NS);
  lmod.set_par(XPar::ktbb, kTbb);
  lmod.set_par(XPar::boost, 0); // direct emission only
  lmod.set_par(XPar::a, spin);
  lmod.set_par(XPar::rin, rin);
  lmod.set_par(XPar::incl, incl);
  VecD flux(n_ener, 0.0);
  XspecSpectrum spec(ener, flux.data(), n_ener);
  lmod.eval_model(spec);
  unsetenv("RELXILL_RRAD_BBODY_TPROFILE");
  return max_shape_deviation(ener, spec.flux, spec_ref, n_ener, emin, emax);
}

// ======= DISK EMISSION: temperature profiles and radial zones (see doc/returnrad_disk_emission.md) ======= //
// (1)  temperature profile and normalization: multicolor disk on a fine radial grid vs. the Xspec models
// (2a) kernel of the direct disk emission vs. the line profiles of the zones (code internal)
// (2b) full relxillBB_NS model with the Novikov-Thorne profile vs. the relativistic Xspec model kerrbb
// (3)  direct and reflected radiation from the same emitted primary (code internal)
// run all of them with bin/tests "[diskemission]"

// temperature profile of the disk emission, and the Xspec model with the same temperature profile
struct DiskProfileCase {
  int tprofile;
  const char *name;
  const char *xspec_model; // also the column name in the reference file
  std::string (*fname_refdata)(double kT); // reference file of the Xspec model
  double max_dev_xspec; // (1) allowed shape deviation of the fine radial grid from Xspec
};

static const DiskProfileCase disk_profile_alpha{
  TPROFILE_ALPHA, "alpha disk (zero torque)", "ezdiskbb",
  get_fname_refdata_ezdiskbb, 0.003
};
static const DiskProfileCase disk_profile_diskbb{
  TPROFILE_DISKBB, "diskbb", "diskbb",
  get_fname_refdata_diskbb, 0.0005
};

// ========== //
// (1) The multicolor disk (spec_multicolor_disk, i.e., the temperature profile and the blackbody spectra of the model,
// without relativistic effects) on a fine radial grid agrees with the Xspec model of the same temperature profile, and
// its normalization changes with kT as for the Xspec model. The outer radius (1e4 Rin) is large enough not to affect
// the band 0.1-10 kT. For the alpha disk, the remaining deviation of ~0.2% (at 10 kT) is due to the radius of the
// maximal temperature, 1.5^(4/5) Rin instead of 49/36 Rin (for which the deviation would be 0.01%).
TEST_CASE (
"disk emission (1): temperature profile on a fine radial grid vs. Xspec diskbb and ezdiskbb"
,
"[rrad-bbody][diskemission]"
)
 {
  const auto profile = GENERATE(disk_profile_alpha, disk_profile_diskbb);

  // grid of the reference files
  constexpr int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  VecD rlo, rhi;
  get_log_radial_grid(1.0, 1e4, 1000, rlo, rhi);

  const std::array<double, 3> kT_values{0.5, 1.0, 2.0};
  std::array<double, 3> norm{};
  for (size_t ii = 0; ii < kT_values.size(); ii++) {
    const double kT = kT_values[ii];
    int status = EXIT_SUCCESS;
    const auto spec = calc_disk_spectrum_rebinned(rlo, rhi, kT, profile.tprofile, ener, n_ener, &status);
    const auto spec_xspec =
        fits_read_refdata_xspec(profile.fname_refdata(kT), profile.xspec_model, ener, n_ener, &status);
    REQUIRE(status == EXIT_SUCCESS);
    const double dev = max_shape_deviation(ener, spec.data(), spec_xspec, n_ener, 0.1 * kT, 10 * kT, &norm[ii]);
    INFO(profile.name << " vs. Xspec " << profile.xspec_model << ", kT = " << kT << ": max. shape deviation "
      << dev << " (allowed: " << profile.max_dev_xspec << ")");
    CHECK(dev < profile.max_dev_xspec);
  }

  // the normalization relative to the Xspec model does not depend on kT
  INFO(profile.name << ": normalization relative to Xspec " << profile.xspec_model << " for kT = 0.5, 2 keV "
    << "(relative to kT = 1 keV): " << norm[0] / norm[1] << ", " << norm[2] / norm[1]);
  CHECK(norm[0] / norm[1] == Catch::Approx(1.0).epsilon(1e-3));
  CHECK(norm[2] / norm[1] == Catch::Approx(1.0).epsilon(1e-3));
}

// ========== //
// (2a) The observed direct disk emission of the model is a single convolution of a blackbody with a kernel:
// the relativistic line profile of the whole disk, with the temperature of each radial point as its line energy
// (calc_relline_profile_radial_line_energy, see get_disk_emission_kernel). With a line energy of 1 keV and the same
// weight for all radial points, this profile is identical to the sum of the line profiles of the zones of the model,
// and a constant line energy E0 shifts it by the factor E0 in energy (without changing the number of photons).
TEST_CASE (
"disk emission (2a): kernel, i.e., line profile with a line energy for each radial point"
,
"[rrad-bbody][diskemission]"
)
 {
  LocalModel lmod(ModelName::relxillBB_NS);
  lmod.set_par(XPar::a, 0.998);
  relParam *rel_param = get_rel_params(lmod.get_model_params());
  xillParam *xill_param = get_xill_params(lmod.get_model_params());
  EnerGrid *egrid = get_relxill_conv_energy_grid();
  const int n_ener = egrid->nbins;

  // line profiles of the zones of the model (as in relxillbb_kernel, also with the xillver table, as the profiles are
  // cached for the model)
  int status = EXIT_SUCCESS;
  xillTable *xill_tab = nullptr;
  get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, &status);
  returningFractions *dat = get_rrad_fractions(rel_param->a, rel_param->rin, rel_param->rout, &status);
  VecD zone_edges(dat->rlo);
  zone_edges.push_back(dat->rhi.back());
  RelSysPar *sys_par = get_system_parameters(rel_param, &status);
  const relline_spec_multizone *profile_zones =
      relbase_profile(egrid->ener, n_ener, rel_param, sys_par, xill_tab, zone_edges, &status);
  REQUIRE(status == EXIT_SUCCESS);

  const int nr = sys_par->nr;
  const VecD weight(nr, 1.0);
  const auto calc_profile = [&](double line_ener) {
    const VecD line_ener_points(nr, line_ener);
    VecD profile(n_ener);
    calc_relline_profile_radial_line_energy(egrid->ener,
                                            n_ener,
                                            sys_par,
                                            line_ener_points.data(),
                                            weight.data(),
                                            profile.data(),
                                            &status);
    return profile;
  };
  const VecD profile_1keV = calc_profile(1.0);
  constexpr double line_ener_shifted = 0.5;
  const VecD profile_shifted = calc_profile(line_ener_shifted);
  REQUIRE(status == EXIT_SUCCESS);

  // identical to the sum of the line profiles of the zones
  double max_abs_dev = 0.0;
  double max_profile = 0.0;
  for (int jj = 0; jj < n_ener; jj++) {
    double sum_zones = 0.0;
    for (int izone = 0; izone < profile_zones->n_zones; izone++) {
      sum_zones += profile_zones->flux[izone][jj];
    }
    max_abs_dev = fmax(max_abs_dev, fabs(profile_1keV[jj] - sum_zones));
    max_profile = fmax(max_profile, sum_zones);
  }
  INFO("line energy 1 keV: max. deviation from the sum of the zone profiles (relative to the maximum): "
    << max_abs_dev / max_profile);
  CHECK(max_abs_dev < 1e-10 * max_profile);

  // a line energy of E0 shifts the profile by the factor E0 (mean of log(E)) and keeps the number of photons
  double sum_1keV = 0.0, sum_shifted = 0.0, mean_log_1keV = 0.0, mean_log_shifted = 0.0;
  for (int jj = 0; jj < n_ener; jj++) {
    const double log_emean = log(0.5 * (egrid->ener[jj] + egrid->ener[jj + 1]));
    sum_1keV += profile_1keV[jj];
    sum_shifted += profile_shifted[jj];
    mean_log_1keV += profile_1keV[jj] * log_emean;
    mean_log_shifted += profile_shifted[jj] * log_emean;
  }
  mean_log_1keV /= sum_1keV;
  mean_log_shifted /= sum_shifted;
  INFO("line energy " << line_ener_shifted << " keV: ratio of the number of photons " << sum_shifted / sum_1keV
    << ", shift of the mean of log(E) " << mean_log_shifted - mean_log_1keV << " (expected "
    << log(line_ener_shifted) << ")");
  CHECK(sum_shifted / sum_1keV == Catch::Approx(1.0).epsilon(1e-3));
  CHECK(mean_log_shifted - mean_log_1keV == Catch::Approx(log(line_ener_shifted)).epsilon(1e-3));

  delete_returningFractions(&dat);
  delete_rel_params(rel_param);
  delete xill_param;
}

// ========== //
// (2b) The direct disk emission (boost=0) of the full relxillBB_NS model with the Novikov-Thorne temperature profile
// (RELXILL_RRAD_BBODY_TPROFILE=3) compared to the relativistic disk model kerrbb of Xspec (eta=0, hd=1, rflag=0,
// lflag=0), with the maximal temperature of the kerrbb disk (1 keV) as kTbb, i.e., without any free parameter. This
// tests the implementation of the disk emission in the model (temperature profile, relativistic effects, and its
// kernel, see (2a)). The spins are on the grid of the kerrbb table, to avoid its linear interpolation in spin.
// The shape deviation is 0.2-1.4% in 0.5-10 kTbb, which is the precision of the kerrbb table (a converged calculation
// on 400 radial zones deviates by the same amount). Lower energies are not compared, as the outer radius of the model
// (1000 Rg) truncates the disk emission there (16% at 0.2 kTbb for a=0). The alpha-disk profile (default) deviates by
// 2.5% or more, i.e., the test is sensitive to the temperature profile.
TEST_CASE (
"disk emission (2b): full relxillBB_NS model with the Novikov-Thorne profile vs. Xspec kerrbb"
,
"[rrad-bbody][diskemission][kerrbb]"
)
 {
  const double spin = GENERATE(0.0, 0.904277, 0.997676); // on the spin grid of kerrbb (no interpolation in spin)
  const double incl = GENERATE(30.0, 60.0);
  constexpr double kTbb = 1.0; // maximal temperature of the disk of the kerrbb reference spectra
  constexpr double rin = -1; // ISCO
  constexpr double max_dev_allowed = 0.02;

  // grid of the reference files
  constexpr int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  int status = EXIT_SUCCESS;
  const auto kerrbb_xspec = fits_read_refdata_xspec(get_fname_refdata_kerrbb(spin, incl),
                                                    "kerrbb",
                                                    ener,
                                                    n_ener,
                                                    &status);
  REQUIRE(status == EXIT_SUCCESS);

  const double emin = 0.5 * kTbb;
  const double emax = 10 * kTbb;
  const double dev_novikov_thorne =
      shape_deviation_relxillbb_ns_disk_emission("3", kTbb, spin, rin, incl, ener, n_ener, kerrbb_xspec, emin, emax);
  const double dev_alpha =
      shape_deviation_relxillbb_ns_disk_emission("1", kTbb, spin, rin, incl, ener, n_ener, kerrbb_xspec, emin, emax);
  INFO("a = " << spin << ", incl = " << incl << ": max. shape deviation from Xspec kerrbb: Novikov-Thorne profile "
    << dev_novikov_thorne << ", alpha-disk profile " << dev_alpha << " (allowed: " << max_dev_allowed << ")");
  CHECK(dev_novikov_thorne < max_dev_allowed);
  CHECK(dev_alpha > max_dev_allowed);
}

// ======= (3) DIRECT AND REFLECTED RADIATION FROM THE SAME EMITTED PRIMARY ======= //
// The reflection is driven by the returning radiation (specRet), which is calculated from the primary emitted by each
// zone (specPri), while the direct emission is calculated with the kernel from the radial points of the line profiles
// (see (2a)). Direct and reflected radiation are therefore only consistent if both are calculated from the same
// emitted primary:
// (3a) specPri of each zone is the emission of the radial points of the line profiles in the zone, and all zones
//      together emit the same as the disk on a fine radial grid
// (3b) specRet is the returning radiation of exactly this primary
// (3c) the primary of the zones, observed with the line profiles of the zones, gives the same flux as the direct
//      emission of the model

// temperature profile (and the value of RELXILL_RRAD_BBODY_TPROFILE) and spin of the tests (3)
struct PrimaryCase {
  const char *env_tprofile;
  int tprofile;
  double spin;
};

#define GENERATE_PRIMARY_CASES GENERATE(PrimaryCase{"1", TPROFILE_ALPHA, 0.0},                \
                                        PrimaryCase{"1", TPROFILE_ALPHA, 0.998},              \
                                        PrimaryCase{"3", TPROFILE_NOVIKOV_THORNE, 0.0},       \
                                        PrimaryCase{"3", TPROFILE_NOVIKOV_THORNE, 0.998})

// emission of the N_FRAD radial points of the line profiles (get_fine_radial_grid, uniform in 1/sqrt(r)) binned into
// the zones of dat, with the area of each point and the split of its cell at the zone boundaries as in
// calc_relline_profile: photon flux per keV and unit area of each zone, zero outside [EMIN_XILLVER, EMAX_XILLVER]
static VecD_2D calc_primary_zones_radial_points(const double *ener,
                                                int n_ener,
                                                const returningFractions *dat,
                                                double rin,
                                                double rout,
                                                double kT,
                                                int tprofile,
                                                int *status) {
  constexpr int nr = N_FRAD;
  VecD re(nr);
  get_fine_radial_grid(rin, rout, re.data(), nr); // descending
  VecD radius(nr); // ascending, as TPROFILE_DISKBB needs the inner edge first
  for (int ii = 0; ii < nr; ii++) {
    radius[ii] = fmax(re[nr - 1 - ii], rin);
  }
  const VecD temperature = get_tprofile(radius, radius, nr, rin, kT, dat->a, tprofile, status);

  VecD zone_edges(dat->rlo);
  zone_edges.push_back(dat->rhi.back());
  VecD_2D spec(dat->nrad, VecD(n_ener, 0.0));
  VecD area(dat->nrad, 0.0);
  VecD spec_point(n_ener);
  for (int ii = 0; ii < nr; ii++) {
    bbody_photon_flux_per_keV(ener, n_ener, spec_point.data(), temperature[nr - 1 - ii], 1.0);
    // cell of the point, split at the zone boundaries in proportion to the area
    const double r_cell_hi = (ii == 0) ? re[ii] : 0.5 * (re[ii - 1] + re[ii]);
    const double r_cell_lo = (ii == nr - 1) ? re[ii] : 0.5 * (re[ii] + re[ii + 1]);
    for (int izone = 0; izone < dat->nrad; izone++) {
      const double r_part_lo = fmax(r_cell_lo, zone_edges[izone]);
      const double r_part_hi = fmin(r_cell_hi, zone_edges[izone + 1]);
      if (r_part_hi > r_part_lo) {
        const double area_part = trapez_integ_single(re.data(), ii, nr) * (
                                   r_part_hi * r_part_hi - r_part_lo * r_part_lo)
                                 / (r_cell_hi * r_cell_hi - r_cell_lo * r_cell_lo);
        area[izone] += area_part;
        for (int jj = 0; jj < n_ener; jj++) {
          spec[izone][jj] += area_part * spec_point[jj];
        }
      }
    }
  }
  for (int ii = 0; ii < dat->nrad; ii++) {
    for (int jj = 0; jj < n_ener; jj++) {
      const bool outside = ener[jj] < EMIN_XILLVER || ener[jj + 1] > EMAX_XILLVER;
      spec[ii][jj] = outside ? 0.0 : spec[ii][jj] / area[ii];
    }
  }
  return spec;
}

// number of photons and energy of a bin-integrated spectrum in the band [emin, emax]
static std::array<double, 2> photons_and_energy(const double *ener,
                                                const double *spec,
                                                int n_ener,
                                                double emin,
                                                double emax) {
  std::array<double, 2> sum{0.0, 0.0};
  for (int jj = 0; jj < n_ener; jj++) {
    if (ener[jj] >= emin && ener[jj + 1] <= emax) {
      sum[0] += spec[jj];
      sum[1] += spec[jj] * 0.5 * (ener[jj] + ener[jj + 1]);
    }
  }
  return sum;
}

// ========== //
// (3a) The primary of each zone (specPri) is the emission of the radial points of the line profiles in the zone
// (identical calculation), and all zones together, weighted with their area, emit the same number of photons and the
// same energy as the disk on a fine radial grid (5000 logarithmic zones).
TEST_CASE (
"disk emission (3a): primary of the zones is the emission of the radial points of the line profiles"
,
"[rrad-bbody][diskemission]"
)
 {
  const auto pcase = GENERATE_PRIMARY_CASES;
  setenv("RELXILL_RRAD_BBODY_TPROFILE", pcase.env_tprofile, 1);
  constexpr double kT = 1.0;
  const double rin = kerr_rms(pcase.spin);
  constexpr double rout = RMAX_RELRET;

  // energy grid of the primary in the model
  int status = EXIT_SUCCESS;
  int n_ener;
  double *ener;
  get_std_bbody_energy_grid(&n_ener, &ener, &status);
  returnSpec2D *return_spec = spec_returnrad_blackbody(ener, n_ener, kT, rin, rout, pcase.spin, &status);
  returningFractions *dat = get_rrad_fractions(pcase.spin, rin, rout, &status);
  const VecD_2D primary_points =
      calc_primary_zones_radial_points(ener, n_ener, dat, rin, rout, kT, pcase.tprofile, &status);
  REQUIRE(status == EXIT_SUCCESS);

  // each zone: the emission of its radial points (bin-integrated)
  double max_dev_zones = 0.0;
  for (int ii = 0; ii < dat->nrad; ii++) {
    double max_spec = 0.0;
    double max_diff = 0.0;
    for (int jj = 0; jj < n_ener; jj++) {
      const double spec_points = primary_points[ii][jj] * (ener[jj + 1] - ener[jj]);
      max_spec = fmax(max_spec, spec_points);
      max_diff = fmax(max_diff, fabs(return_spec->specPri[ii][jj] - spec_points));
    }
    max_dev_zones = fmax(max_dev_zones, max_diff / max_spec);
  }

  // all zones together vs. the disk on a fine radial grid
  VecD rlo_fine;
  VecD rhi_fine;
  get_log_radial_grid(rin, rout, 5000, rlo_fine, rhi_fine);
  const VecD temp_fine =
      get_tprofile(rlo_fine, rhi_fine, static_cast<int>(rlo_fine.size()), rin, kT, pcase.spin, pcase.tprofile, &status);
  VecD spec_fine(n_ener, 0.0);
  VecD spec_zone(n_ener);
  for (size_t kk = 0; kk < rlo_fine.size(); kk++) {
    bbody_photon_flux_per_keV(ener, n_ener, spec_zone.data(), temp_fine[kk], 1.0);
    for (int jj = 0; jj < n_ener; jj++) {
      spec_fine[jj] += M_PI * (rhi_fine[kk] * rhi_fine[kk] - rlo_fine[kk] * rlo_fine[kk]) * spec_zone[jj]
          * (ener[jj + 1] - ener[jj]);
    }
  }
  VecD spec_all_zones(n_ener, 0.0);
  for (int ii = 0; ii < dat->nrad; ii++) {
    for (int jj = 0; jj < n_ener; jj++) {
      spec_all_zones[jj] += M_PI * (dat->rhi[ii] * dat->rhi[ii] - dat->rlo[ii] * dat->rlo[ii])
          * return_spec->specPri[ii][jj];
    }
  }
  const auto sum_zones = photons_and_energy(ener, spec_all_zones.data(), n_ener, EMIN_XILLVER, EMAX_XILLVER);
  const auto sum_fine = photons_and_energy(ener, spec_fine.data(), n_ener, EMIN_XILLVER, EMAX_XILLVER);
  REQUIRE(status == EXIT_SUCCESS);

  INFO("tprofile " << pcase.tprofile << ", a = " << pcase.spin << ": max. deviation of specPri of a zone from its "
    << "radial points " << max_dev_zones << "; all zones vs. fine radial grid: photons "
    << sum_zones[0] / sum_fine[0] - 1 << ", energy " << sum_zones[1] / sum_fine[1] - 1);
  CHECK(max_dev_zones < 1e-10);
  CHECK(fabs(sum_zones[0] / sum_fine[0] - 1) < 1e-3);
  CHECK(fabs(sum_zones[1] / sum_fine[1] - 1) < 1e-3);

  unsetenv("RELXILL_RRAD_BBODY_TPROFILE");
  delete_returnSpec2D(&return_spec);
  delete_returningFractions(&dat);
}

// ========== //
// (3b) The returning radiation of each zone (specRet) is the returning radiation of exactly the primary of (3a).
TEST_CASE (
"disk emission (3b): returning radiation of the zones is calculated from this primary"
,
"[rrad-bbody][diskemission]"
)
 {
  const auto pcase = GENERATE_PRIMARY_CASES;
  setenv("RELXILL_RRAD_BBODY_TPROFILE", pcase.env_tprofile, 1);
  constexpr double kT = 1.0;
  const double rin = kerr_rms(pcase.spin);
  constexpr double rout = RMAX_RELRET;

  int status = EXIT_SUCCESS;
  int n_ener;
  double *ener;
  get_std_bbody_energy_grid(&n_ener, &ener, &status);
  returnSpec2D *return_spec = spec_returnrad_blackbody(ener, n_ener, kT, rin, rout, pcase.spin, &status);
  returningFractions *dat = get_rrad_fractions(pcase.spin, rin, rout, &status);
  const int norder = get_from_env_rrad_bbody_order();
  const VecD_2D primary_points =
      calc_primary_zones_radial_points(ener, n_ener, dat, rin, rout, kT, pcase.tprofile, &status);
  const VecD_2D returnrad_points = get_returnrad_specs(ener, n_ener, dat, primary_points, norder, &status);

  REQUIRE(status == EXIT_SUCCESS);

  double max_dev_points = 0.0;
  for (int ii = 0; ii < dat->nrad; ii++) {
    const double energy_model = photons_and_energy(ener, return_spec->specRet[ii], n_ener, 0.0, 1e10)[1];
    const double energy_points = photons_and_energy(ener, returnrad_points[ii].data(), n_ener, 0.0, 1e10)[1];
    double max_spec = 0.0;
    double max_diff = 0.0;
    for (int jj = 0; jj < n_ener; jj++) {
      max_spec = fmax(max_spec, returnrad_points[ii][jj]);
      max_diff = fmax(max_diff, fabs(return_spec->specRet[ii][jj] - returnrad_points[ii][jj]));
    }
    max_dev_points = fmax(max_dev_points, fmax(max_diff / max_spec, fabs(energy_model / energy_points - 1)));
  }
  INFO("tprofile " << pcase.tprofile << ", a = " << pcase.spin << ": max. deviation of specRet from the returning "
    << "radiation of the primary of the radial points " << max_dev_points);
  CHECK(max_dev_points < 1e-10);

  unsetenv("RELXILL_RRAD_BBODY_TPROFILE");
  delete_returnSpec2D(&return_spec);
  delete_returningFractions(&dat);
}

// ========== //
// (3c) The direct disk emission of the full relxillBB_NS model (boost=0, i.e., the kernel of (2a)) and the primary of
// the zones (specPri, from which the reflection is calculated, see (3a) and (3b)), observed with the line profiles of
// the zones, have the same number of photons and the same energy flux (band 0.1-1000 keV, as for the normalization of
// the reflection). They only differ by the variation of the temperature within a zone combined with the relativistic
// effects, which is not taken into account for the zones (0.02-0.11% in photons and 0.05-0.17% in energy). This
// requires that the line profiles of the zones cover exactly the area of the zones, as the returning radiation
// (calc_relline_profile splits the cells of the radial points at the zone boundaries; assigning each cell completely
// to one zone instead changes the area of a zone by up to 4%, which results in a deviation of 0.4-1.2% here).
TEST_CASE (
"disk emission (3c): direct emission of the model vs. the primary of the zones observed"
,
"[rrad-bbody][diskemission]"
)
 {
  const auto pcase = GENERATE_PRIMARY_CASES;
  setenv("RELXILL_RRAD_BBODY_TPROFILE", pcase.env_tprofile, 1);
  constexpr double kTbb = 1.0;

  LocalModel lmod(ModelName::relxillBB_NS);
  lmod.set_par(XPar::a, pcase.spin);
  lmod.set_par(XPar::incl, 30.0);
  lmod.set_par(XPar::rin, -1);
  lmod.set_par(XPar::ktbb, kTbb);
  lmod.set_par(XPar::boost, 0); // direct emission only

  // direct emission of the model, on the energy grid of the convolution
  EnerGrid *egrid = get_relxill_conv_energy_grid();
  const int n_ener = egrid->nbins;
  VecD flux_model(n_ener, 0.0);
  XspecSpectrum spec_model(egrid->ener, flux_model.data(), n_ener);
  lmod.eval_model(spec_model);

  // primary of the zones convolved with the line profiles of the zones (as in relxillbb_kernel)
  int status = EXIT_SUCCESS;
  relParam *rel_param = get_rel_params(lmod.get_model_params());
  xillParam *xill_param = get_xill_params(lmod.get_model_params());
  xillTable *xill_tab = nullptr;
  get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, &status);
  returnSpec2D *return_spec =
      spec_returnrad_blackbody(egrid->ener, n_ener, kTbb, rel_param->rin, rel_param->rout, rel_param->a, &status);
  RelSysPar *sys_par = get_system_parameters(rel_param, &status);
  const relline_spec_multizone *profiles = relbase_profile(egrid->ener,
                                                           n_ener,
                                                           rel_param,
                                                           sys_par,
                                                           xill_tab,
                                                           get_radial_grid_from_returntab(return_spec),
                                                           &status);
  specCache *spec_cache = init_global_specCache(&status);
  VecD flux_zones(n_ener, 0.0);
  VecD flux_conv(n_ener);
  for (int ii = 0; ii < return_spec->nrad; ii++) {
    convolveSpectrumFFTNormalized(egrid->ener,
                                  return_spec->specPri[ii],
                                  profiles->flux[ii],
                                  flux_conv.data(),
                                  n_ener,
                                  1,
                                  1,
                                  ii,
                                  spec_cache,
                                  &status);
    for (int jj = 0; jj < n_ener; jj++) {
      flux_zones[jj] += flux_conv[jj];
    }
  }
  REQUIRE(status == EXIT_SUCCESS);

  const auto sum_model = photons_and_energy(egrid->ener, flux_model.data(), n_ener, 0.1, 1000.0);
  const auto sum_zones = photons_and_energy(egrid->ener, flux_zones.data(), n_ener, 0.1, 1000.0);
  INFO("tprofile " << pcase.tprofile << ", a = " << pcase.spin << ": primary of the zones observed relative to the "
    << "direct emission of the model (0.1-1000 keV): photons " << sum_zones[0] / sum_model[0] - 1 << ", energy "
    << sum_zones[1] / sum_model[1] - 1);
  CHECK(fabs(sum_zones[0] / sum_model[0] - 1) < 3e-3);
  CHECK(fabs(sum_zones[1] / sum_model[1] - 1) < 3e-3);

  unsetenv("RELXILL_RRAD_BBODY_TPROFILE");
  delete_returnSpec2D(&return_spec);
  delete_rel_params(rel_param);
  delete xill_param;
}
