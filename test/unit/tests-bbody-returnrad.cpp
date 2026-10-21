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
#include "tests-bbody-returnrad.h"

extern "C" {
#include "relutility.h"
}

#if !defined(RELXILL_SOURCE_DIR)
#warning: RELXILL_SOURCE_DIR is not defined
#endif

#define PREC 1e-6
#define PREC_REFDATA 1e-3
#define REFDATA_DIR "test/unit/refdata/"  // relativ to RELXILL_SOURCE_DIR (set when compiling)


// Function definitions
static double *integSpecArea(double **spec, int n_ener, double spin, const returnSpec2D *returnSpec, int *status);

static void fits_write_spec(const char *fname, const double *ener, double *spec, int n_ener, int *status) {
  CHECK_STATUS_VOID(*status);

  char filename[300];
  if (fname[0] != '!') {
    sprintf(filename, "!%s", fname);
  } else {
    sprintf(filename, "%s", fname);
  }

  // open the fits file
  fitsfile *fptr;
  if (fits_create_file(&fptr, filename, status)) {
    CHECK_STATUS_VOID(*status);
    printf("  *** error : creating file %s failed\n", fname);
    return;
  }

  char tf1[3] = "1D";
  char f1[10] = "bin_lo";
  char f2[10] = "bin_hi";
  char f3[10] = "flux";

  const int tfields = 3;
  char *ttype[] = {f1, f2, f3};
  char *tform[] = {tf1, tf1, tf1};

  if (fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "spec", status)) {
    CHECK_STATUS_VOID(*status);
  }

  double bin_lo[n_ener];
  double bin_hi[n_ener];

  for (int ii = 0; ii < n_ener; ii++) {
    bin_lo[ii] = ener[ii];
    bin_hi[ii] = ener[ii + 1];
  }

  int firstrow = 1; /* first row in table to write   */
  int firstelem = 1; /* first element in row  */

  fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, n_ener, bin_lo, status);
  fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, n_ener, bin_hi, status);
  fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, n_ener, spec, status);

  relxill_check_fits_error(status);

  if (fptr != nullptr) { fits_close_file(fptr, status); }
}


static refSpecData *new_refSpecData(int nener, int *status) {
  CHECK_STATUS_RET(*status, nullptr);

  auto *spec = new refSpecData;
  CHECK_MALLOC_RET_STATUS(spec, status, nullptr);

  spec->nener = nener;

  spec->elo = new double[nener];
  spec->ehi = new double[nener];
  spec->flux = new double[nener];

  return spec;
}

static void free_refSpecData(refSpecData *spec) {
  delete[] spec->elo;
  delete[] spec->ehi;
  delete[] spec->flux;
  delete spec;
}

static refSpecData *fits_read_refdata_spec(char *fname, int *status) {
  // open the reference data
  fitsfile *fptr = nullptr;
  (fits_open_table(&fptr, fname, READONLY, status));
  relxill_check_fits_error(status);
  CHECK_STATUS_RET(*status, nullptr);

  long nener;

  // get the number of rows
  if (fits_get_num_rows(fptr, &nener, status)) return nullptr;

  auto spec = new_refSpecData((int) nener, status);
  CHECK_MALLOC_RET_STATUS(spec, status, nullptr);

  int anynul = 0;
  double nullptrval = 0.0;
  auto nelem = (LONGLONG) nener;
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, nelem, &nullptrval, spec->elo, &anynul, status);
  fits_read_col(fptr, TDOUBLE, 2, 1, 1, nelem, &nullptrval, spec->ehi, &anynul, status);
  fits_read_col(fptr, TDOUBLE, 3, 1, 1, nelem, &nullptrval, spec->flux, &anynul, status);

  if (*status != EXIT_SUCCESS) {
    printf(" *** error *** initializing of the RETURN RADIATION table %s failed \n", fname);
    free_refSpecData(spec);
  }

  if (fptr != nullptr) {
    fits_close_file(fptr, status);
  }

  return spec;
}

static int compare_refSpecDat(const refSpecData *spec, const double *ener0, const double *flux0, int nener0) {
  int status = EXIT_SUCCESS;

  if (spec->nener != nener0) {
    status = EXIT_FAILURE;
    printf(" *** error : spectrum has %i bins, but expected to have %i bins\n", spec->nener, nener0);
    return EXIT_FAILURE;
  }

  for (int ii = 0; ii < nener0; ii++) {
    if ((spec->flux[ii] / flux0[ii]) > PREC_REFDATA && (fabs(spec->flux[ii] - flux0[ii]) > PREC_REFDATA)) {
      status = EXIT_FAILURE;
      printf(" *** error : spectrum flux in bin [%.3f, %.3f] is %e, but expected to be %e \n",
             ener0[ii],
             ener0[ii + 1],
             flux0[ii],
             spec->flux[ii]);
    }
  }

  return status;
}

static int compareWithRefdata(const char *fname_ref,
                              const double *ener0,
                              const double *flux0,
                              int nener0,
                              int *status) {
  PRINT_RELXILL_TEST_MSG(fname_ref);

  char fullname_ref[1000];
  sprintf(fullname_ref, "%s/%s/%s", RELXILL_SOURCE_DIR, REFDATA_DIR, fname_ref);
  refSpecData *spec = fits_read_refdata_spec(fullname_ref, status);

  if (spec == nullptr) {
    return -1;
  }

  int status_comp = compare_refSpecDat(spec, ener0, flux0, nener0);
  free_refSpecData(spec);
  return status_comp;
}


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
// "DiskBB" tests below (created with test/unit/refdata/create_refvalues_diskbb.sl)
static std::string get_fname_refdata_diskbb(double Tin) {
  if (fabs(Tin - 1.0) < 1e-6) {
    return "refvalues_diskbb.fits";
  }
  char fname[100];
  snprintf(fname, sizeof(fname), "refvalues_diskbb_Tin%.1f.fits", Tin);
  return fname;
}

// reference file with the ezdiskbb model of Xspec (column "ezdiskbb") for the given Tmax, on the energy grid of the
// "ezdiskbb" tests below (created with test/unit/refdata/create_refvalues_diskbb.sl)
static std::string get_fname_refdata_ezdiskbb(double Tmax) {
  char fname[100];
  snprintf(fname, sizeof(fname), "refvalues_ezdiskbb_Tmax%.1f.fits", Tmax);
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
// spec/spec_ref normalized to its mean in the band
static double max_shape_deviation(const double *ener,
                                  const double *spec,
                                  const VecD &spec_ref,
                                  int n_ener,
                                  double emin,
                                  double emax) {
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
  return max_dev;
}

// multicolor disk spectrum (spec_diskbb or spec_ezdiskbb) on the grid ener: as in the model, it is calculated on the
// fine energy grid of the convolution and then rebinned (evaluating the blackbody at the bin centres of the coarse
// grid directly would add a deviation of ~3%)
static VecD calc_disk_spectrum_rebinned(void (*spec_disk)(const double *, double *, int, double, double, int *),
                                        double kT,
                                        double spin,
                                        const double *ener,
                                        int n_ener,
                                        int *status) {
  const EnerGrid *conv_grid = get_relxill_conv_energy_grid();
  VecD spec_conv(conv_grid->nbins);
  spec_disk(conv_grid->ener, spec_conv.data(), conv_grid->nbins, kT, spin, status);
  VecD spec(n_ener, 0.0);
  _rebin_spectrum(ener, spec.data(), n_ener, conv_grid->ener, spec_conv.data(), conv_grid->nbins);
  return spec;
}

// max. shape deviation (see max_shape_deviation) of the direct disk emission (boost=0) of the full relxillBB_NS model
// from spec_ref in the band [emin, emax], for the temperature profile tprofile (value of RELXILL_RRAD_BBODY_TPROFILE)
static double shape_deviation_relxillbb_ns_disk_emission(const char *tprofile,
                                                         double kTbb,
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

// ========== //
TEST_CASE (

" DiskBB Spectrum "
,
"[bbret]"
)
{

  const double Tin = GENERATE(0.5, 1.0, 2.0);

  int status = EXIT_SUCCESS;
  setenv("DEBUG_RELXILL","1",1);

  //create an energy grid (grid of the reference files)
  constexpr int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  double spin = 0.998; // simply determines the radial grid here

  auto photar = calc_disk_spectrum_rebinned(spec_diskbb, Tin, spin, ener, n_ener, &status);
  REQUIRE(status == EXIT_SUCCESS);

  if (fabs(Tin - 1.0) < 1e-6) {
    fits_write_spec("!testrr-spec-diskbb.fits", ener, photar.data(), n_ener, &status);
  }

  // the shape agrees with the diskbb model of Xspec; the deviation of ~1.6% is due to the 50 radial zones (mainly
  // the coarse outer zones, which affect the lowest energies), see doc/returnrad_disk_emission.md
  const auto diskbb_xspec = fits_read_refdata_xspec(get_fname_refdata_diskbb(Tin), "diskbb", ener, n_ener, &status);
  REQUIRE(status == EXIT_SUCCESS);
  const double dev = max_shape_deviation(ener, photar.data(), diskbb_xspec, n_ener, 0.1 * Tin, 10 * Tin);
  INFO("Tin = " << Tin << ": max. shape deviation from Xspec diskbb = " << dev);
  REQUIRE(dev < 0.02);

  unsetenv("DEBUG_RELXILL");
}

// ========== //
TEST_CASE(" ezdiskbb Spectrum ", "[bbret]") {
  const double Tmax = GENERATE(0.5, 1.0, 2.0);
  int status = EXIT_SUCCESS;

  // grid of the reference files
  constexpr int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  const double spin = 0.998; // simply determines the radial grid here

  const auto photar = calc_disk_spectrum_rebinned(spec_ezdiskbb, Tmax, spin, ener, n_ener, &status);
  REQUIRE(status == EXIT_SUCCESS);

  // the shape agrees with the ezdiskbb model of Xspec (zero-torque profile with Tmax=kTbb, the default profile of the
  // relxillBB models); the deviation of ~4% is due to the 50 radial zones (mainly the coarse outer zones, which affect
  // the lowest energies; a continuous disk agrees to 0.3%), see doc/returnrad_disk_emission.md
  const auto ezdiskbb_xspec =
      fits_read_refdata_xspec(get_fname_refdata_ezdiskbb(Tmax), "ezdiskbb", ener, n_ener, &status);
  REQUIRE(status == EXIT_SUCCESS);
  const double dev = max_shape_deviation(ener, photar.data(), ezdiskbb_xspec, n_ener, 0.1 * Tmax, 10 * Tmax);
  INFO("Tmax = " << Tmax << ": max. shape deviation from Xspec ezdiskbb = " << dev);
  REQUIRE(dev < 0.05);
}

// ========== //
// Direct disk emission (boost=0) of the full relxillBB_NS model (temperature profile, blackbody spectra of the
// zones, relativistic convolution) compared to the diskbb model of Xspec:
// - diskbb temperature profile (RELXILL_RRAD_BBODY_TPROFILE=2): the shape agrees to ~8%; the remaining difference
//   is relativistic, and for Rin=10 r_isco and incl=50deg gravitational redshift and Doppler shift partly cancel
//   (other values of Rin=-8...-12 and incl=40...50deg give deviations of 7-20%), see doc/returnrad_disk_emission.md
// - zero-torque temperature profile (default): the shape differs by ~30%, i.e., the test is sensitive to the
//   temperature profile
TEST_CASE(" relxillBB_NS direct disk emission vs. Xspec diskbb for the diskbb and zero-torque temperature profiles ",
          "[bbret][diskbb]") {
  const double Tin = GENERATE(0.5, 1.0, 2.0);
  constexpr double max_dev_allowed = 0.1;
  constexpr double rin = -10;
  constexpr double incl = 50;

  // grid of the reference files
  constexpr int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  int status = EXIT_SUCCESS;
  const auto diskbb_xspec = fits_read_refdata_xspec(get_fname_refdata_diskbb(Tin), "diskbb", ener, n_ener, &status);
  REQUIRE(status == EXIT_SUCCESS);

  const double emin = 0.15 * Tin;
  const double emax = 10 * Tin;
  const double dev_diskbb =
      shape_deviation_relxillbb_ns_disk_emission("2", Tin, rin, incl, ener, n_ener, diskbb_xspec, emin, emax);
  const double dev_zero_torque =
      shape_deviation_relxillbb_ns_disk_emission("1", Tin, rin, incl, ener, n_ener, diskbb_xspec, emin, emax);
  INFO("Tin = " << Tin << ": max. shape deviation from Xspec diskbb: diskbb profile " << dev_diskbb
       << ", zero-torque profile " << dev_zero_torque << " (allowed: " << max_dev_allowed << ")");
  CHECK(dev_diskbb < max_dev_allowed);
  CHECK(dev_zero_torque > max_dev_allowed);
}

// ========== //
// Direct disk emission (boost=0) of the full relxillBB_NS model compared to the ezdiskbb model of Xspec:
// - zero-torque temperature profile (default): the shape agrees to ~7% in 0.3-10 Tmax for Rin=10 r_isco and
//   incl=40deg. Lower energies are not compared: for this Rin, the disk only extends over a factor of 80 in radius
//   (Rout=1000Rg), which truncates the emission of the outer disk (relatively hotter than for diskbb), see
//   doc/returnrad_disk_emission.md
// - diskbb temperature profile: the shape differs by ~30%, i.e., the test is sensitive to the temperature profile
TEST_CASE(" relxillBB_NS direct disk emission vs. Xspec ezdiskbb for the zero-torque and diskbb temperature profiles ",
          "[bbret][ezdiskbb]") {
  const double Tmax = GENERATE(0.5, 1.0, 2.0);
  constexpr double max_dev_allowed = 0.1;
  constexpr double rin = -10;
  constexpr double incl = 40;

  // grid of the reference files
  constexpr int n_ener = 50;
  double ener[n_ener + 1];
  get_log_grid(ener, n_ener + 1, 0.01, 20.0);

  int status = EXIT_SUCCESS;
  const auto ezdiskbb_xspec =
      fits_read_refdata_xspec(get_fname_refdata_ezdiskbb(Tmax), "ezdiskbb", ener, n_ener, &status);
  REQUIRE(status == EXIT_SUCCESS);

  const double emin = 0.3 * Tmax;
  const double emax = 10 * Tmax;
  const double dev_zero_torque =
      shape_deviation_relxillbb_ns_disk_emission("1", Tmax, rin, incl, ener, n_ener, ezdiskbb_xspec, emin, emax);
  const double dev_diskbb =
      shape_deviation_relxillbb_ns_disk_emission("2", Tmax, rin, incl, ener, n_ener, ezdiskbb_xspec, emin, emax);
  INFO("Tmax = " << Tmax << ": max. shape deviation from Xspec ezdiskbb: zero-torque profile " << dev_zero_torque
    << ", diskbb profile " << dev_diskbb << " (allowed: " << max_dev_allowed << ")");
  CHECK(dev_zero_torque < max_dev_allowed);
  CHECK(dev_diskbb > max_dev_allowed);
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

  // refdata currently missing
  // REQUIRE( compareWithRefdata("refvalues_rrbbody_ret.fits", ener, photar_areaInteg, n_ener, &status) == EXIT_SUCCESS);
  // REQUIRE( compareWithRefdata("refvalues_rrbbody_pri.fits", ener, photar0_areaInteg, n_ener, &status ) == EXIT_SUCCESS);

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
