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
}

#if !defined(RELXILL_SOURCE_DIR)
#warning: RELXILL_SOURCE_DIR is not defined
#endif

#define PREC 1e-6
#define REFDATA_DIR "test/unit/refdata/"  // relativ to RELXILL_SOURCE_DIR (set when compiling)


// Function definitions
static double *integSpecArea(double **spec, int n_ener, double spin, const returnSpec2D *returnSpec, int *status);

static double *integSpecArea(double **spec, int n_ener, double spin, const returnSpec2D *returnSpec, int *status) {
  auto *areaIntegSpec = (double *) malloc(sizeof(double) * n_ener);
  CHECK_MALLOC_RET_STATUS(areaIntegSpec, status, areaIntegSpec)

  for (int jj = 0; jj < n_ener; jj++) {
    areaIntegSpec[jj] = 0.0;
  }
  for (int ii = 0; ii < returnSpec->nrad; ii++) {
    double prop_area = calc_proper_area_ring(returnSpec->rlo[ii], returnSpec->rhi[ii], spin);
    for (int jj = 0; jj < n_ener; jj++) {
      areaIntegSpec[jj] += spec[ii][jj] * prop_area;
    }
  }
  return areaIntegSpec;
}


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
  snprintf(fname, sizeof(fname), "refvalues_kerrbb_a%.3f_i%.0f.fits", spin, incl);
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

static int checkTemperatureProfile(double *temperature, int n) {
  if (temperature[0] < 0.0) {
    printf("\n *** error: lowest bin of temperature profile is %f  and NOT >=0 as expected \n",
           temperature[0]);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

void fitsAddKeyRadialZone(char *fname, const returnSpec2D *returnSpec, int izone, int *status) {
  fitsfile *fptr = nullptr;
  fits_open_table(&fptr, fname, READWRITE, status);

  double rlo = returnSpec->rlo[izone];
  double rhi = returnSpec->rhi[izone];
  fits_update_key(fptr, TDOUBLE, "rlo", &(rlo), "radial zone grid", status);
  fits_update_key(fptr, TDOUBLE, "rhi", &(rhi), "radial zone grid", status);
  fits_update_key(fptr, TINT, "izone", &izone, "radial zone grid", status);

  fits_close_file(fptr, status);


  //relxill_check_fits_error(status);
}


//  ======= TEST FUNCTIONS ========   //

// ========== //
TEST_CASE (

" Check Temperature Profile"
,
"[bbret]"
)
{

  int status = EXIT_SUCCESS;

  std::vector<double> spinArray = { 0.85,0.851,0.86,0.87,0.89,0.9 };

  double Rout = 1000;

  returningFractions* dat = nullptr;
  for (double spin : spinArray){
    double Rin = kerr_rms(spin);
    dat = get_rrad_fractions(spin, Rin, Rout, &status);

   VecD temperature = get_temperature_profile_disk_zones(dat, Rin, 1.0, &status);

    REQUIRE( checkTemperatureProfile(temperature.data(), dat->nrad) == EXIT_SUCCESS);

    delete_returningFractions(&dat);
  }
}

// ======= DISK EMISSION: temperature profiles and radial zones (see doc/returnrad_disk_emission.md) ======= //
// (1)  temperature profile and normalization: multicolor disk on a fine radial grid vs. the Xspec models
// (2a) radial zones of the model vs. a fine radial grid (code internal)
// (2b) full relxillBB_NS model with the Novikov-Thorne profile vs. the relativistic Xspec model kerrbb
// run all of them with bin/tests "[diskemission]"

// temperature profile of the disk emission, and the Xspec model with the same temperature profile
struct DiskProfileCase {
  int tprofile;
  const char *name;
  const char *xspec_model; // also the column name in the reference file
  std::string (*fname_refdata)(double kT); // reference file of the Xspec model
  double max_dev_xspec; // (1) allowed shape deviation of the fine radial grid from Xspec
  double max_dev_zones; // (2a) allowed shape deviation of the model zones from the fine grid
};

static const DiskProfileCase disk_profile_alpha{
  TPROFILE_ALPHA, "alpha disk (zero torque)", "ezdiskbb",
  get_fname_refdata_ezdiskbb, 0.003, 0.04
};
static const DiskProfileCase disk_profile_diskbb{
  TPROFILE_DISKBB, "diskbb", "diskbb",
  get_fname_refdata_diskbb, 0.0005, 0.02
};

// ========== //
// (1) The multicolor disk (spec_multicolor_disk, i.e., the temperature profile and the blackbody spectra of the model,
// without relativistic effects) on a fine radial grid agrees with the Xspec model of the same temperature profile, and
// its normalization changes with kT as for the Xspec model. The outer radius (1e4 Rin) is large enough not to affect
// the band 0.1-10 kT. For the alpha disk, the remaining deviation of ~0.2% (at 10 kT) is due to the radius of the
// maximal temperature, 1.5^(4/5) Rin instead of 49/36 Rin (for which the deviation would be 0.01%).
TEST_CASE(" disk emission (1): temperature profile on a fine radial grid vs. Xspec diskbb and ezdiskbb ",
          "[bbret][diskemission]") {
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
// (2a) The disk emission on the radial zones of the model (from the returning radiation table, 50 zones from r_isco to
// 1000 Rg) compared to a fine radial grid over the same radii. This is a code internal comparison, as (1) shows that
// the fine radial grid is correct. The deviation (alpha disk ~3.5%, diskbb ~1.6% in 0.1-10 kT, mainly at the lowest
// energies from the coarse outer zones) is the precision of the disk emission of the model without relativistic
// effects. It does not depend on kT, as the shape of the spectrum only depends on E/kT.
TEST_CASE(" disk emission (2a): radial zones of the model vs. a fine radial grid ", "[bbret][diskemission]") {
  const auto profile = GENERATE(disk_profile_alpha, disk_profile_diskbb);
  constexpr double kT = 1.0;
  constexpr double spin = 0.998;

  constexpr int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  int status = EXIT_SUCCESS;
  returningFractions *dat = get_rrad_fractions(spin, kerr_rms(spin), RMAX_RELRET, &status);
  REQUIRE(status == EXIT_SUCCESS);
  const auto spec_zones = calc_disk_spectrum_rebinned(dat->rlo, dat->rhi, kT, profile.tprofile, ener, n_ener, &status);

  VecD rlo, rhi;
  get_log_radial_grid(dat->rlo.front(), dat->rhi.back(), 500, rlo, rhi);
  const auto spec_fine = calc_disk_spectrum_rebinned(rlo, rhi, kT, profile.tprofile, ener, n_ener, &status);
  REQUIRE(status == EXIT_SUCCESS);

  const double dev = max_shape_deviation(ener, spec_zones.data(), spec_fine, n_ener, 0.1 * kT, 10 * kT);
  INFO(profile.name << ": " << dat->nrad << " radial zones of the model vs. a fine radial grid: max. shape deviation "
    << dev << " (allowed: " << profile.max_dev_zones << ")");
  CHECK(dev < profile.max_dev_zones);

  delete_returningFractions(&dat);
}

// ========== //
// (2b) The direct disk emission (boost=0) of the full relxillBB_NS model with the Novikov-Thorne temperature profile
// (RELXILL_RRAD_BBODY_TPROFILE=3) compared to the relativistic disk model kerrbb of Xspec (eta=0, hd=1, rflag=0,
// lflag=0), with the maximal temperature of the kerrbb disk (1 keV) as kTbb, i.e., without any free parameter. This
// tests the implementation of the disk emission in the model (temperature profile, radial zones, relativistic
// effects). The shape deviation is 0.1-1% in 0.5-10 kTbb. Lower energies are not compared, as the outer radius of the
// model (1000 Rg) truncates the disk emission there (15% at 0.2 kTbb for a=0). The alpha-disk profile (default)
// deviates by more than 4%, i.e., the test is sensitive to the temperature profile.
TEST_CASE(" disk emission (2b): full relxillBB_NS model with the Novikov-Thorne profile vs. Xspec kerrbb ",
          "[bbret][diskemission][kerrbb]") {
  const double spin = GENERATE(0.0, 0.9, 0.998);
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

// ========== //
TEST_CASE (

" Return Radiation Black Body Spectrum"
,
"[bbret]"
)
{

  int status = EXIT_SUCCESS;

  /* create an energy grid */
  int n_ener = 200;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 50.0);

  /* call the relline model */
  double photar[n_ener];
  double photar0[n_ener];

  double Tin = 1.0;
  double spin = 0.998; // simply determines the radial grid here

  double Rin = kerr_rms(spin);
  double Rout = RMAX_RELRET;

  returnSpec2D *returnSpec = spec_returnrad_blackbody(ener, n_ener, Tin, Rin, Rout, spin, &status);

  assert(returnSpec != nullptr);

  double* photar_areaInteg = integSpecArea(returnSpec->specRet, n_ener, spin, returnSpec, &status);
  double* photar0_areaInteg = integSpecArea(returnSpec->specPri, n_ener, spin, returnSpec, &status);

  fits_rr_write_2Dspec("!testrr-rframe-rr-bbody.fits", returnSpec->specRet, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, &status);
  fits_rr_write_2Dspec("!testrr-rframe-prim-bbody.fits", returnSpec->specPri, ener, n_ener,
                       returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, &status);

  free(photar_areaInteg);
  free(photar0_areaInteg);

  delete_returnSpec2D(&returnSpec);

}

TEST_CASE (

" Check loading of the xillverNS table for relxillBB_NS"
,
"[bbret]"
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




TEST_CASE (

"Check if Second Evaluation returns the identical results"
,
"[bbret]"
)
{
  LocalModel lmod(ModelName::relxillBB);

  DefaultSpec default_spec1{};
  XspecSpectrum spec1 = default_spec1.get_xspec_spectrum();
  lmod.eval_model(spec1);

  DefaultSpec default_spec0{};
  XspecSpectrum spec0 = default_spec0.get_xspec_spectrum();
  lmod.eval_model(spec0);

  auto sum0 = calcSum(spec0.flux,spec0.num_flux_bins());
  REQUIRE( std::isnormal (sum0));
  auto sum1 = calcSum(spec1.flux,spec0.num_flux_bins());
    REQUIRE(fabs(sum0 - sum1) < PREC);
}

static void test_boost_parameter(LocalModel &lmod) {
  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();

  lmod.set_par(XPar::boost, 1);
  lmod.eval_model(spec);
  double sum_total = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_total > 1e-8);

  lmod.set_par(XPar::boost, -1);
  lmod.eval_model(spec);
  double sum_refl = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_refl > 1e-8);

  lmod.set_par(XPar::boost, 0);
  lmod.eval_model(spec);
  double sum_direct = calcSum(spec.flux, spec.num_flux_bins());
  REQUIRE(sum_direct > 1e-8);


  // check that boost parameter works: refl+direct has to give total flux
  INFO(" total: " << sum_total << " direct: " << sum_direct << " refl: " << sum_refl
    << " (diff: " << fabs(sum_total - sum_direct - sum_refl) / sum_total << ")");
  REQUIRE(fabs(sum_total - sum_direct - sum_refl) / sum_total < 1e-4);
}


TEST_CASE(" relxillBB boost parameter: direct + reflected = total ", "[relxillBB]") {
  LocalModel lmod(ModelName::relxillBB);
  test_boost_parameter(lmod);
}
