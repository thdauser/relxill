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

#include <algorithm>
#include <iostream>
#include <numeric>
#include "Relreturn_BlackBody.h"
#include "ModelDefinition.h"
#include "Rebin.h"
#include "Relbase.h"
#include "Relphysics.h"
#include "Relreturn_Datastruct.h"

#include "writeOutfiles.h"

extern "C" {
#include "xilltable_reading.h"
#include "common.h"
#include "relutility.h"
}

#define LIM_GFAC_RR_BBODY 0.001 // difference between gmin and gmax, above which the energy shift is taken into account
const std::string fbase_rrad_bbody = "!__relxillOutput-rrad-bbody-";  // filename base for output files
const std::string fbase_rrad_bbody_dat = "__relxillOutput-rrad-bbody-"; // filename base for output files


// Function definitions
void bin_integrate_flux(int nener, const double* ener, VecD_2D& spec);

// Program Code


void fits_rr_write_2Dspec(const char* fname, double** spec_arr, const double* ener, int nener,
                          const VecD& rlo, const VecD& rhi, int nrad, returningFractions* dat, int* status)
{
    CHECK_STATUS_VOID(*status);

    // open the fits file
    fitsfile* fptr;
    if (fits_create_file(&fptr, fname, status))
    {
        relxill_check_fits_error(status);
        printf("   creating file %s failed\n", fname);
        CHECK_STATUS_VOID(*status);
    }


    int n1 = nrad;
    int n2 = nener;
    char dim1[10];
    char dim2[10];
    sprintf(dim1, "%iD", 1);
    sprintf(dim2, "%iD", n2);

    if (dat != nullptr)
    {
        constexpr int tfields = 6;
        char *ttype[] = {
            const_cast<char *>("rlo"), const_cast<char *>("rhi"), const_cast<char *>("ener"),
            const_cast<char *>("spec"), const_cast<char *>("fraci"), const_cast<char *>("fret")
        };
        char *tform[] = {dim1, dim1, dim2, dim2, const_cast<char *>("50D"), dim1};
        fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);
    } else {
        constexpr int tfields = 4;
        char *ttype[] = {
            const_cast<char *>("rlo"), const_cast<char *>("rhi"),
            const_cast<char *>("ener"), const_cast<char *>("spec")};
        char* tform[] = {dim1, dim1, dim2, dim2};
        fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, nullptr, "2Dspec", status);
    }

    relxill_check_fits_error(status);
    CHECK_STATUS_VOID(*status);

    int firstrow = 1; /* first row in table to write   */
    int firstelem = 1; /* first element in a row */

    // we use "const_cast" as with legacy cfitsio this is the only way (and assume those values are never changed)
    fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, n1, const_cast<double*>(rlo.data()), status);
    fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, n1, const_cast<double*>(rhi.data()), status);
    for (int ii = 0; ii < n1; ii++)
    {
        fits_write_col(fptr, TDOUBLE, 3, firstrow + ii, firstelem, n2, const_cast<double*>(ener), status);
        fits_write_col(fptr, TDOUBLE, 4, firstrow + ii, firstelem, n2, spec_arr[ii], status);
    }

    if (dat != nullptr)
    {
        for (int ii = 0; ii < n1; ii++)
        {
            fits_write_col(fptr, TDOUBLE, 5, firstrow + ii, firstelem, n1, dat->tf_r[ii].data(), status);
        }
        fits_write_col(fptr, TDOUBLE, 6, firstrow, firstelem, n1, dat->tabData->f_ret, status);
        if (fits_write_key(fptr, TDOUBLE, "SPIN", &(dat->a), nullptr, status))
        {
        }
    }


    relxill_check_fits_error(status);

    if (fptr != nullptr) { fits_close_file(fptr, status); }
}


void fits_rr_write_2Dspec(const char* fname, const std::vector<std::vector<double>>& spec_arr,
                          const double* ener, int nener,
                          const VecD& rlo, const VecD& rhi, int nrad, returningFractions* dat, int* status)
{
    // Convert to double** temporarily
    std::vector<double*> ptrs(spec_arr.size());
    for (size_t i = 0; i < spec_arr.size(); i++)
    {
        ptrs[i] = const_cast<double*>(spec_arr[i].data());
    }

    // Call the original implementation
    fits_rr_write_2Dspec(fname, ptrs.data(), ener, nener, rlo, rhi, nrad, dat, status);
}

// temperature profile of the disk emission, selected with the environment variable RELXILL_RRAD_BBODY_TPROFILE
static int get_tprofile_type_from_env() {
    switch (get_from_env_rrad_bbody_tprofile()) {
        case 2: return TPROFILE_DISKBB;
        case 3: return TPROFILE_NOVIKOV_THORNE;
        default: return TPROFILE_ALPHA; // this is for a value of 1 as well
    }
}

// temperature of each radial zone (evaluated at its mean radius), only used for the diagnostic output
VecD get_temperature_profile_disk_zones(const returningFractions *dat, double Rin, double Tin, int *status) {
    return get_tprofile(dat->rlo, dat->rhi, dat->nrad, Rin, Tin, dat->a, get_tprofile_type_from_env(), status);
}

// temperature of the profile tprofile (see get_tprofile) at the radial points re, which are descending as for the
// line profiles (get_fine_radial_grid); the innermost point is rin up to rounding errors
static VecD get_temperature_radial_points(const double *re,
                                          int nr,
                                          double rin,
                                          double kT,
                                          double spin,
                                          int tprofile,
                                          int *status) {
    assert(re[0] >= re[nr - 1]);
    // ascending radii, as TPROFILE_DISKBB needs the inner edge as first radius
    VecD radius(nr);
    for (int ii = 0; ii < nr; ii++) {
        radius[ii] = fmax(re[nr - 1 - ii], rin);
    }
    const VecD temperature = get_tprofile(radius, radius, nr, rin, kT, spin, tprofile, status);
    return {temperature.rbegin(), temperature.rend()}; // same (descending) order as re
}


/**
 * @brief Emitted primary (direct) disk emission of each zone, from the radial points of the line profiles.
 *
 * The disk is integrated over the same N_FRAD radial points as the relativistic line profiles (get_fine_radial_grid
 * from rin to rout, uniform in 1/sqrt(r)), with the same area of each point and the same split of its cell at the zone
 * boundaries as in calc_relline_profile, and with the same temperature of each point as the direct disk emission of
 * the model (get_disk_emission_kernel). The primary of a zone is the mean of the blackbodies of its points, weighted
 * with the area of their cells in the zone. As the returning radiation is calculated from this primary, the reflected
 * and the direct emission of the model are based on the same emitted primary (see the unit tests "disk emission
 * (3)").
 *
 * @param ener Energy bin boundaries [keV] (size nener+1).
 * @param nener Number of energy bins.
 * @param rlo Inner radii of the zones [r_g].
 * @param rhi Outer radii of the zones [r_g].
 * @param rin Inner edge of the disk [r_g].
 * @param rout Outer edge of the disk [r_g].
 * @param kT Temperature of the profile [keV] (see get_tprofile).
 * @param spin Spin (only used for TPROFILE_NOVIKOV_THORNE).
 * @param tprofile Temperature profile (TPROFILE_ALPHA, TPROFILE_DISKBB, or TPROFILE_NOVIKOV_THORNE).
 * @param status Status flag.
 *
 * @return Photon flux per keV and unit area of each zone [nzones][nener] (not bin-integrated), set to zero outside
 *         [EMIN_XILLVER, EMAX_XILLVER].
 */
static VecD_2D get_disk_emission_zones(const double *ener,
                                       int nener,
                                       const VecD &rlo,
                                       const VecD &rhi,
                                       double rin,
                                       double rout,
                                       double kT,
                                       double spin,
                                       int tprofile,
                                       int *status) {
    CHECK_STATUS_RET(*status, VecD_2D());
    const int nzones = static_cast<int>(rlo.size());
    assert(nzones > 0 && rhi.size() == rlo.size());

    VecD rgrid(rlo);
    rgrid.push_back(rhi.back());

    VecD re(N_FRAD);
    get_fine_radial_grid(rin, rout, re.data(), N_FRAD); // descending, as for the line profiles
    const VecD temperature = get_temperature_radial_points(re.data(), N_FRAD, rin, kT, spin, tprofile, status);
    CHECK_STATUS_RET(*status, VecD_2D());

    // sum of the blackbodies of the points of each zone, weighted with their area; the blackbody is only calculated up
    // to E = 60 kT (above, it is below exp(-60) of its maximum)
    VecD_2D spec_zones(nzones, VecD(nener, 0.0));
    VecD area_zones(nzones, 0.0);
    VecD spec_point(nener);
    for (int ii = 0; ii < N_FRAD; ii++) {
        const int nener_bb = std::min(
            nener,
            static_cast<int>(std::upper_bound(ener, ener + nener, 60 * temperature[ii]) - ener));
        if (temperature[ii] > 0) {
            // zero at the inner edge for zero torque
            bbody_photon_flux_per_keV(ener, nener_bb, spec_point.data(), temperature[ii], 1.0);
        }

        // as in calc_relline_profile, the cell of the point (between the midpoints to its neighbors) is split at the
        // boundaries of the zones in proportion to the area (radii outside of rgrid belong to the first or last zone)
        const double r_cell_hi = (ii == 0) ? re[ii] : 0.5 * (re[ii - 1] + re[ii]);
        const double r_cell_lo = (ii == N_FRAD - 1) ? re[ii] : 0.5 * (re[ii] + re[ii + 1]);
        const double area_cell = trapez_integ_single(re.data(), ii, N_FRAD);
        const int izone_lo = binary_search(rgrid.data(), nzones + 1, r_cell_lo);
        const int izone_hi = binary_search(rgrid.data(), nzones + 1, r_cell_hi);
        for (int izone = izone_lo; izone <= izone_hi; izone++) {
            const double r_part_lo = (izone == 0) ? r_cell_lo : fmax(r_cell_lo, rgrid[izone]);
            const double r_part_hi = (izone == nzones - 1) ? r_cell_hi : fmin(r_cell_hi, rgrid[izone + 1]);
            if (r_part_hi <= r_part_lo) {
                continue;
            }
            const double area = area_cell * (r_part_hi * r_part_hi - r_part_lo * r_part_lo)
                                / (r_cell_hi * r_cell_hi - r_cell_lo * r_cell_lo);
            area_zones[izone] += area;
            if (temperature[ii] > 0) {
                for (int jj = 0; jj < nener_bb; jj++) {
                    spec_zones[izone][jj] += area * spec_point[jj];
                }
            }
        }
    }

    // average over the area of each zone, and zero outside [EMIN_XILLVER, EMAX_XILLVER]
    for (int ii = 0; ii < nzones; ii++) {
        assert(area_zones[ii] > 0); // the zones and the cells of the points both cover the disk from rin to rout
        for (int jj = 0; jj < nener; jj++) {
            const bool outside = ener[jj] < EMIN_XILLVER || ener[jj + 1] > EMAX_XILLVER;
            spec_zones[ii][jj] = outside ? 0.0 : spec_zones[ii][jj] / area_zones[ii];
        }
    }
    return spec_zones;
}

/**
 * @brief Direct (primary) disk emission of the given radial zones, calculated on a fine radial grid.
 *
 * Each zone [rlo, rhi] is divided into logarithmic sub-zones with r_hi/r_lo <= RRAD_BBODY_MAX_RATIO_SUBZONES. Each
 * sub-zone emits a blackbody with the temperature of the profile tprofile at its mean radius, and the spectrum of the
 * zone is the average of these blackbodies, weighted with the (flat-space) area of the sub-zones. This gives the
 * correct spectrum also for coarse zones, for which a single temperature is not a good approximation (e.g., the outer
 * zones of the returning-radiation table, see doc/returnrad_disk_emission.md).
 *
 * @param ener Energy bin boundaries [keV] (size nener+1).
 * @param nener Number of energy bins.
 * @param rlo Inner radii of the zones [r_g].
 * @param rhi Outer radii of the zones [r_g].
 * @param Rin Inner edge of the disk [r_g].
 * @param kT Temperature of the profile [keV]: maximal temperature (TPROFILE_ALPHA, TPROFILE_NOVIKOV_THORNE), or
 *           temperature at Rin (TPROFILE_DISKBB).
 * @param spin Spin (only used for TPROFILE_NOVIKOV_THORNE).
 * @param tprofile Temperature profile (TPROFILE_ALPHA, TPROFILE_DISKBB, or TPROFILE_NOVIKOV_THORNE).
 * @param status Status flag.
 *
 * @return Bin-integrated spectra per unit area [nzones][nener], set to zero outside [EMIN_XILLVER, EMAX_XILLVER].
 */
VecD_2D get_disk_emission_zones_binned(const double *ener,
                                       int nener,
                                       const VecD &rlo,
                                       const VecD &rhi,
                                       double Rin,
                                       double kT,
                                       double spin,
                                       int tprofile,
                                       int *status) {
    CHECK_STATUS_RET(*status, VecD_2D());
    const int nzones = static_cast<int>(rlo.size());
    assert(nzones > 0 && rhi.size() == rlo.size());

    // logarithmic sub-zones of each zone
    VecD rlo_sub;
    VecD rhi_sub;
    std::vector<int> izone_sub;
    for (int ii = 0; ii < nzones; ii++) {
        const int nsub =
                std::max(1,
                         static_cast<int>(ceil(log(rhi[ii] / rlo[ii]) / log(RRAD_BBODY_MAX_RATIO_SUBZONES) - 1e-10)));
        for (int kk = 0; kk < nsub; kk++) {
            rlo_sub.push_back(rlo[ii] * pow(rhi[ii] / rlo[ii], static_cast<double>(kk) / nsub));
            rhi_sub.push_back(rlo[ii] * pow(rhi[ii] / rlo[ii], static_cast<double>(kk + 1) / nsub));
            izone_sub.push_back(ii);
        }
    }
    const int nsub_total = static_cast<int>(rlo_sub.size());

    VecD temperature_sub = get_tprofile(rlo_sub, rhi_sub, nsub_total, Rin, kT, spin, tprofile, status);
    CHECK_STATUS_RET(*status, VecD_2D());

    // sum of the blackbody spectra of the sub-zones of each zone, weighted with their area; the blackbody is only
    // calculated up to E = 60 kT (above, it is below exp(-60) of its maximum)
    VecD_2D spec_zones(nzones, VecD(nener, 0.0));
    VecD area_zones(nzones, 0.0);
    VecD spec_sub(nener);
    for (int kk = 0; kk < nsub_total; kk++) {
        const int nener_bb = static_cast<int>(
            std::upper_bound(ener, ener + nener, 60 * temperature_sub[kk]) - ener); // bins below 60 kT
        bbody_photon_flux_per_keV(ener, std::min(nener_bb, nener), spec_sub.data(), temperature_sub[kk], 1.0);

        const double area = M_PI * (rhi_sub[kk] * rhi_sub[kk] - rlo_sub[kk] * rlo_sub[kk]);
        area_zones[izone_sub[kk]] += area;
        for (int jj = 0; jj < std::min(nener_bb, nener); jj++) {
            spec_zones[izone_sub[kk]][jj] += area * spec_sub[jj];
        }
    }

    // average over the area of each zone, bin-integrated, and zero outside [EMIN_XILLVER, EMAX_XILLVER]
    for (int ii = 0; ii < nzones; ii++) {
        for (int jj = 0; jj < nener; jj++) {
            if (ener[jj] < EMIN_XILLVER || ener[jj + 1] > EMAX_XILLVER) {
                spec_zones[ii][jj] = 0.0;
            } else {
                spec_zones[ii][jj] *= (ener[jj + 1] - ener[jj]) / area_zones[ii];
            }
        }
    }
    return spec_zones;
}


/**
 * @brief Relativistic kernel of the direct disk emission, for a maximal temperature of 1 keV (kTbb = 1).
 *
 * Per unit area, the bin-integrated blackbody photon spectrum of the temperature T is T^3 times the one of 1 keV,
 * shifted in energy by the factor T. With the temperature profile T(r) = kTbb t(r), the direct emission of the disk is
 * therefore the convolution of the blackbody of kTbb with a kernel, which is the relativistic line profile of each
 * radial point at the line energy t(r), weighted with t(r)^3 (calc_relline_profile_radial_line_energy). The kernel is
 * integrated over the same N_FRAD radial points as the line profiles of the zones, but with the temperature of each
 * point instead of a single temperature per zone (see doc/returnrad_disk_emission.md).
 *
 * The kernel does not depend on kTbb. It is cached and only re-calculated if the relativistic parameters or the
 * temperature profile change.
 *
 * @param ener Energy grid of the convolution [keV] (size n_ener+1).
 * @param n_ener Number of energy bins.
 * @param rel_param Relativistic parameters.
 * @param sysPar System parameters of rel_param (radial points, transfer functions).
 * @param tprofile Temperature profile (TPROFILE_ALPHA, TPROFILE_DISKBB, or TPROFILE_NOVIKOV_THORNE).
 * @param status Status flag.
 *
 * @return Kernel [photons/bin] (size n_ener), to be convolved with the bin-integrated blackbody photon spectrum
 *         (bbody_photon_flux_per_keV) per unit area.
 */
static const VecD &get_disk_emission_kernel(const double *ener,
                                            int n_ener,
                                            const relParam *rel_param,
                                            RelSysPar *sysPar,
                                            int tprofile,
                                            int *status) {
    static relParam *cached_rel_param = nullptr;
    static int cached_tprofile = -1;
    static VecD cached_kernel;

    if (did_rel_param_change(cached_rel_param, rel_param) == 0 && tprofile == cached_tprofile
        && cached_kernel.size() == static_cast<size_t>(n_ener)) {
        return cached_kernel;
    }
    cached_tprofile = -1; // invalid until the new kernel is calculated

    // line energy: temperature of each radial point for kTbb = 1 (as for the primary of the zones,
    // get_disk_emission_zones)
    const int nr = sysPar->nr;
    const VecD line_ener =
            get_temperature_radial_points(sysPar->re, nr, rel_param->rin, 1.0, rel_param->a, tprofile, status);
    CHECK_STATUS_RET(*status, cached_kernel);

    VecD weight(nr);
    for (int ii = 0; ii < nr; ii++) {
        weight[ii] = pow(line_ener[ii], 3);
    }

    cached_kernel.assign(n_ener, 0.0);
    calc_relline_profile_radial_line_energy(ener,
                                            n_ener,
                                            sysPar,
                                            line_ener.data(),
                                            weight.data(),
                                            cached_kernel.data(),
                                            status);
    CHECK_STATUS_RET(*status, cached_kernel);

    set_cached_rel_param(rel_param, &cached_rel_param, status);
    cached_tprofile = tprofile;
    return cached_kernel;
}


/**
 * @brief Observed direct disk emission of the relxillBB models [photons/bin], on the energy grid of the convolution.
 *
 * Convolution of the blackbody of the maximal temperature kTbb (bin-integrated photon spectrum per unit area) with the
 * kernel of the disk emission (get_disk_emission_kernel). The FFT convolution is cyclic, and the kernel shifts the
 * energies down by up to a factor ~100 (outer disk). The blackbody is therefore set to zero in the bins that the
 * kernel can shift outside the energy grid, as they would re-appear at the other end of the grid. This only changes
 * the spectrum observed below ~0.05 keV (~0.1 keV for the diskbb profile).
 */
static VecD calc_disk_emission_observed(double *ener,
                                        int n_ener,
                                        double kTbb,
                                        const relParam *rel_param,
                                        RelSysPar *sysPar,
                                        specCache *spec_cache,
                                        int *status) {
    VecD spec_obs(n_ener, 0.0);
    const VecD &kernel =
            get_disk_emission_kernel(ener, n_ener, rel_param, sysPar, get_tprofile_type_from_env(), status);
    CHECK_STATUS_RET(*status, spec_obs);

    // maximal shift of the kernel to lower and higher energies (in bins, 1 keV is no shift as in fftw_conv_spectrum)
    const auto is_positive = [](double val) { return val > 0; };
    const int ibin_1keV = binary_search(ener, n_ener + 1, 1.0);
    const int nbins_shift_down =
            ibin_1keV - static_cast<int>(std::find_if(kernel.begin(), kernel.end(), is_positive) - kernel.begin());
    const int nbins_shift_up =
            static_cast<int>(kernel.rend() - std::find_if(kernel.rbegin(), kernel.rend(), is_positive)) - 1 - ibin_1keV;

    VecD spec_bbody(n_ener);
    bbody_photon_flux_per_keV(ener, n_ener, spec_bbody.data(), kTbb, 1.0);
    for (int jj = 0; jj < n_ener; jj++) {
        if (jj < nbins_shift_down || jj >= n_ener - nbins_shift_up) {
            spec_bbody[jj] = 0.0;
        } else {
            spec_bbody[jj] *= ener[jj + 1] - ener[jj]; // bin-integrated
        }
    }

    // the FFT arrays of zone 0 are only used as buffer (the FFTs are always re-calculated in relxillbb_kernel)
    convolveSpectrumFFTNormalized(ener,
                                  spec_bbody.data(),
                                  kernel.data(),
                                  spec_obs.data(),
                                  n_ener,
                                  1,
                                  1,
                                  0,
                                  spec_cache,
                                  status);

    // convolveSpectrumFFTNormalized normalizes the number of photons in the band [EMIN_XILLVER, EMAX_XILLVER], but the
    // cool outer disk is observed with up to ~1% of its photons below EMIN_XILLVER: normalize it on the whole energy
    // grid instead, where the convolution conserves the number of photons (no bins are wrapped around)
    const double sum_obs = std::accumulate(spec_obs.begin(), spec_obs.end(), 0.0);
    const double norm_fac = std::accumulate(kernel.begin(), kernel.end(), 0.0)
                            * std::accumulate(spec_bbody.begin(), spec_bbody.end(), 0.0) / sum_obs;
    for (double &val: spec_obs) {
        val *= norm_fac;
    }
    return spec_obs;
}


/**
 * @brief Calculates a 2D blackbody return radiation spectrum for relativistic disk models
 *
 * This function computes the returning radiation spectrum from a relativistic accretion disk
 * using blackbody emission. The spectrum is calculated across multiple radial zones and
 * energy bins, taking into account relativistic effects near a spinning black hole.
 *
 * @param ener Pointer to array of energy values [keV] (input energy grid)
 * @param nener Number of energy bins in the energy array
 * @param Tin Temperature at the inner disk radius [K]
 * @param Rin Inner disk radius [gravitational radii, r_g = GM/c²]
 * @param Rout Outer disk radius [gravitational radii]
 * @param spin Dimensionless black hole spin parameter (a* = a/M, range: -1 to 1)
 * @param status Pointer to status variable (0 = success, non-zero = error)
 *
 * @return Pointer to returnSpec2D structure containing:
 *         - 2D spectral array [n_radial_zones][n_energy_bins]
 *         - Radial grid information
 *         - Energy grid information
 *         Returns NULL on error
 *
 * @note The temperature profile is selected with RELXILL_RRAD_BBODY_TPROFILE (default: alpha disk with zero torque
 * at the inner edge, see doc/returnrad_disk_emission.md)
 * @note The primary (direct) disk emission of each zone (specPri) is calculated from the radial points of the line
 * profiles, binned into the zones (get_disk_emission_zones), and the returning radiation (specRet) is calculated from
 * exactly this primary. The observed direct emission of the relxillBB models is calculated from the same radial
 * points and temperatures, but for the whole disk at once (calc_disk_emission_observed).
 * @note Memory for the returned structure must be freed by the caller
 * @note The spin parameter follows the convention where positive values indicate prograde rotation
 *
 * @see returnSpec2D for the structure definition
 * @see spec_multicolor_disk for the same disk emission without relativistic effects (used in the unit tests)
 */
returnSpec2D* spec_returnrad_blackbody(double* ener, int nener,
                                       double Tin,
                                       double Rin,
                                       double Rout,
                                       double spin,
                                       int *status) {
    CHECK_STATUS_RET(*status, NULL);

    // 1 - get the fractions from the table (plus interpolation to current parameters)
    returningFractions *dat = get_rrad_fractions(spin, Rin, Rout, status);


    // 2 - temperature of each zone at its mean radius (only for the diagnostic output)
    VecD temperature = get_temperature_profile_disk_zones(dat, Rin, Tin, status);

    // 3 - primary (direct) disk emission of each zone, from the radial points of the line profiles in the zone, on the
    //     energy grid of the returning radiation
    int nener_bbody;
    double *ener_bbody;
    get_std_bbody_energy_grid(&nener_bbody, &ener_bbody, status);
    VecD_2D spec_prim_zones_bbody = get_disk_emission_zones(ener_bbody,
                                                            nener_bbody,
                                                            dat->rlo,
                                                            dat->rhi,
                                                            Rin,
                                                            Rout,
                                                            Tin,
                                                            dat->a,
                                                            get_tprofile_type_from_env(),
                                                            status);

    // 4 - returning radiation of exactly this primary, and both on the output grid (bin integrated spectra)
    int num_rrad_order = get_from_env_rrad_bbody_order();
    VecD_2D spec_rr_zones = get_returnrad_specs(ener, nener, dat, spec_prim_zones_bbody, num_rrad_order, status);
    if (*status != EXIT_SUCCESS) {
        delete_returningFractions(&dat);
        return nullptr;
    }

    VecD_2D spec_prim_zones(dat->nrad, VecD(nener));
    for (int ii = 0; ii < dat->nrad; ii++) {
        rebin_mean_flux(ener, spec_prim_zones[ii].data(), nener, ener_bbody, spec_prim_zones_bbody[ii].data(),
                        nener_bbody, status);
    }
    bin_integrate_flux(nener, ener, spec_prim_zones);

    if (is_debug_run())
    {
        for (int ii=0; ii<dat->nrad; ii++)
        {
            assert_double_array_isnormal(spec_rr_zones[ii].data(), nener);
            assert_double_array_isnormal(spec_prim_zones[ii].data(), nener);
        }
    }

    returnSpec2D *returnSpec = get_returnrad_output_structure(dat, spec_rr_zones, spec_prim_zones, ener, nener, status);

    if (shouldOutfilesBeWritten()) {
        std::string fname = "!undefined.fits";
        fname = fbase_rrad_bbody + "rframe-bbodyEmit.fits";
        fits_rr_write_2Dspec(fname.c_str(), returnSpec->specPri, ener, nener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);

        fname = fbase_rrad_bbody + "rframe-bbodyRet.fits";
        fits_rr_write_2Dspec(fname.c_str(), returnSpec->specRet, ener, nener,
                             returnSpec->rlo, returnSpec->rhi, returnSpec->nrad, nullptr, status);


        fname = fbase_rrad_bbody_dat + "rframe-bbodyTemp.dat";
        FILE *fp = fopen(fname.c_str(), "w+");
        fprintf(fp, "#zone \t rlo \t rhi \t T_col \n");
        for (int ii = 0; ii < returnSpec->nrad; ii++) {
            fprintf(fp,
                    "%i \t %e \t %e \t %e \n",
                    ii,
                    returnSpec->rlo[ii],
                    returnSpec->rhi[
                        +ii],
                    temperature[ii]);
        }
        fclose_errormsg(fp, fname.c_str());
    }

    delete_returningFractions(&dat);

    return returnSpec;
}


/*
 * Helper Routines for Black Body returning radiation calculations
 */


static void calc_rr_bbspec_gzone(double *ener,
                                 int nener,
                                 double *spec_ret,
                                 const VecD &spec_pri,
                                 const double *gfac,
                                 int ng,
                                 const double* frac_g)
{
    // loop over all gfac and calculate the spectrum for each
    VecD spec_ret_g(nener, 0.0);
    setArrayToZero(spec_ret, nener);
    for (int kk = 0; kk < ng; kk++)
    {
        // we have to shift 1/gfac, as gfac=E_i/E_e and do not shift the spectrum, but shift the energy grid "backwards"
        // N_i(E_i) = N_e(E_e) = N_e(E_i / gfac)
        rebin_mean_flux_shifted(ener, nener, spec_ret_g, spec_pri, 1./gfac[kk]);
        for (int jj = 0; jj < nener; jj++)
        {
            spec_ret[jj] += spec_ret_g[jj] * frac_g[kk];
        }
    }
}


/**
 * @brief Compute the returning radiation blackbody spectrum for a single radial zone.
 *
 * Sums contributions from all radial zones that return photons to zone `irad`,
 * applying the appropriate energy shift (g-factor) and returning fractions.
 * The result is a spectrum in ph/s/keV/cm² (not bin-integrated).
 *
 * @param[out] spec     Output spectrum array (size nener), zeroed and filled on return.
 * @param      ener     Energy grid boundaries (size nener+1) [keV].
 * @param      nener    Number of energy bins.
 * @param      irad     Index of the receiving radial zone.
 * @param      spec_pri Primary blackbody spectra per zone (ph/s/keV/cm²).
 * @param      dat      Returning fractions data (g-factors, transfer fractions, f_ret).
 */
static void calc_rr_bbspec_ring(VecD& spec,
                                double* ener,
                                int nener,
                                int irad,
                                const VecD_2D &spec_pri,
                                const returningFractions* dat)
{
    VecD spec_r(nener, 0.0);
    assert(spec.size() == nener);
    spec.assign(spec.size(), 0.0);

    for (int ii = 0; ii < dat->nrad; ii++)
    {
        auto gfac = get_gfac_grid(dat->tabData->gmin[irad][ii], dat->tabData->gmax[irad][ii], dat->tabData->ng);

        // only calculate it for a significant redshift, otherwise return BBODY without gshift
        if (fabs(dat->tabData->gmax[irad][ii] - dat->tabData->gmin[irad][ii]) > LIM_GFAC_RR_BBODY)
        {
            calc_rr_bbspec_gzone(ener, nener, spec_r.data(), spec_pri[ii], gfac.data(), dat->tabData->ng, dat->tabData->frac_g[irad][ii]);
        }
        else
        {
            // we have to shift 1/gfac, as gfac=E_i/E_e and do not shift the spectrum, but shift the energy grid "backwards"
            // N_i(E_i) = N_e(E_e) = N_e(E_i / gfac)
            rebin_mean_flux_shifted(ener, nener, spec_r, spec_pri[ii], 1./dat->tabData->gmin[irad][ii]);
        }

        // apply the correct fractions and add it to the zone output spectrum
        for (int jj = 0; jj < nener; jj++)
        {
            spec[jj] += spec_r[jj] * dat->tabData->f_ret[ii] * dat->tf_r[irad][ii];
        }
    }
    assert_double_array_isnormal(spec.data(), nener);
}


/**
 * @brief Compute returning radiation spectra from the given primary spectra for all zones.
 *
 * Loops over all radial zones and calls `calc_rr_bbspec_ring` to apply the
 * returning fractions and energy shifts. The result is not bin-integrated.
 *
 * @param      dat      Returning fractions data (g-factors, transfer fractions, f_ret).
 * @param      nener    Number of energy bins.
 * @param      ener     Energy grid boundaries (size nener+1) [keV].
 * @param      spec_pri Primary spectra per zone (ph/s/keV/cm²).
 *
 * @return 2D spectrum [nrad][nener] in ph/s/keV/cm² (not bin-integrated).
 */
static VecD_2D calc_returnrad_from_input_spectra(const returningFractions *dat,
                                                 int nener,
                                                 double *ener,
                                                 const VecD_2D &spec_pri) {
    VecD_2D spec_zones(dat->nrad, VecD(nener));
    for (int ii = 0; ii < dat->nrad; ii++) {
        calc_rr_bbspec_ring(spec_zones[ii], ener, nener, ii, spec_pri, dat); // return ph / cm²/s/keV ( not bin integ.)
    }
    return spec_zones;
}


/**
 * @brief Compute bin-integrated returning radiation spectra for all radial zones.
 *
 * Computes the returning radiation of the given primary spectra of the zones (on the internal bbody energy grid,
 * get_std_bbody_energy_grid) via `calc_returnrad_from_input_spectra`, then rebins onto the output grid and
 * bin-integrates.
 *
 * @param      ener_inp    Output energy grid boundaries (size nener_inp+1) [keV].
 * @param      nener_inp   Number of output energy bins.
 * @param      dat         Returning fractions data (g-factors, transfer fractions, f_ret).
 * @param      spec_pri    Primary spectrum of each zone on the bbody energy grid (ph/s/keV/cm², not bin-integrated),
 *                         see get_disk_emission_zones.
 * @param      norder      Number of returning radiation calculations (combined orders are returned!)
 * @param      status      Error status flag.
 *
 * @return 2D spectrum [nrad][nener_inp] in cts/bin, or empty on error.
 */
VecD_2D get_returnrad_specs(double *ener_inp,
                            int nener_inp,
                            const returningFractions *dat,
                            const VecD_2D &spec_pri,
                            int norder,
                            int *status)
{
    CHECK_STATUS_RET(*status, VecD_2D());
    assert(!dat->tf_r[0].empty());

    int nener;
    double* ener;
    get_std_bbody_energy_grid(&nener, &ener, status);
    assert(spec_pri.size() == static_cast<size_t>(dat->nrad) && spec_pri[0].size() == static_cast<size_t>(nener));

    std::vector<VecD_2D> spec_out(norder);
    VecD_2D spec_sum_out(dat->nrad, VecD(nener));
    for (int iorder = 0; iorder < norder; iorder++) {
        spec_out[iorder] = calc_returnrad_from_input_spectra(dat, nener, ener, iorder == 0 ? spec_pri : spec_out[iorder - 1]);
        if (iorder == 0) {
            spec_sum_out = spec_out[iorder];
        } else {
            add_2d(spec_sum_out, spec_out[iorder]);
        }
    }

    // match the spectra to the output grid
    VecD_2D spec_zones_rebinned(dat->nrad, VecD(nener_inp));
    for (int ii = 0; ii < dat->nrad; ii++)
    {
        rebin_mean_flux(ener_inp, spec_zones_rebinned[ii].data(), nener_inp, ener, spec_sum_out[ii].data(), nener, status);
        assert_double_array_isnormal(spec_zones_rebinned[ii].data(),nener_inp);
    }
    bin_integrate_flux(nener_inp, ener_inp, spec_zones_rebinned);

    return spec_zones_rebinned;
}




/**
 * @brief Converts spectral values to counts per energy bin for multiple radial zones
 *
 * This function adjusts the spectral flux values for a given 2D spectra by integrating over
 * energy bins. For each radial zone, the spectral values are multiplied by the energy
 * bin width to produce counts per bin.
 *
 * @param nener Number of energy bins in the input energy array
 * @param ener Pointer to an array containing the boundaries of the energy bins [keV]
 *             (array size must be nener + 1)
 * @param spec 2D vector of spectral values, organized as spec[nrad][nener], containing
 *             the flux values to be converted (e.g., photons/cm²/s/keV).
 */
void bin_integrate_flux(int nener, const double *ener, VecD_2D &spec) {
    const int nrad_spec = static_cast<int>(spec.size());
    for (int ii = 0; ii < nrad_spec; ii++) {
        for (int jj = 0; jj < nener; jj++) {
            spec[ii][jj] *= (ener[jj + 1] - ener[jj]); // now make it cts/bin
        }
    }
}

/**
 * multicolor disk spectrum [photons/bin, arbitrary normalization] for the given radial zones [rlo, rhi]
 *
 * The spectrum of each zone is calculated with get_disk_emission_zones_binned (i.e., on a fine radial grid within each
 * zone), and is weighted with its flat-space ring area (i.e., without relativistic effects):
 * - TPROFILE_ALPHA: alpha disk with zero torque at the inner edge rlo[0], with the maximal temperature kT
 * - TPROFILE_DISKBB: diskbb profile, with the temperature kT at the inner edge rlo[0]
 *
 * It uses the same temperature profiles and blackbody spectra as the disk emission of the relxillBB models, and is
 * used in the unit tests to compare them with the Xspec models ezdiskbb and diskbb (see
 * doc/returnrad_disk_emission.md).
 */
void spec_multicolor_disk(const double *ener,
                          double *spec,
                          int n,
                          const VecD &rlo,
                          const VecD &rhi,
                          double kT,
                          int tprofile,
                          int *status) {
    CHECK_STATUS_VOID(*status);
    const int nrad = static_cast<int>(rlo.size());
    assert(nrad > 0 && rhi.size() == rlo.size());
    assert(tprofile != TPROFILE_NOVIKOV_THORNE); // relativistic profile (needs the spin), not for this function

    VecD_2D spec_arr = get_disk_emission_zones_binned(ener, n, rlo, rhi, rlo[0], kT, 0.0, tprofile, status);
    CHECK_STATUS_VOID(*status);

    for (int jj = 0; jj < n; jj++) {
        spec[jj] = 0.0;
        for (int ii = 0; ii < nrad; ii++) {
            spec[jj] += spec_arr[ii][jj] * M_PI * (rhi[ii] * rhi[ii] - rlo[ii] * rlo[ii]);
        }
    }
}


VecD get_radial_grid_from_returntab(const returnSpec2D *spec) {
    assert(spec != nullptr);
    assert(spec->nrad > 0);
    assert(!spec->rlo.empty());
    assert(!spec->rhi.empty());

    VecD rgrid;
    rgrid.reserve(spec->nrad + 1); // we use n+1 grid points

    // Add lower bounds
    for (int ii = 0; ii < spec->nrad; ii++)
    {
        rgrid.push_back(spec->rlo[ii]);
    }

    // Add the upper bound of the last zone
    rgrid.push_back(spec->rhi[spec->nrad - 1]);

    return rgrid;
}

static double get_xillver_norm_factor_from_primary_spectrum(const double *spec_in,
                                                            const double *ener,
                                                            int n_ener,
                                                            int *status) {
    CHECK_STATUS_RET(*status, 0.0);

    EnerGrid* egrid = get_coarse_xillver_energrid();

    auto xillverInputSpec = new double[egrid->nbins];
    CHECK_MALLOC_RET_STATUS(xillverInputSpec, status, 0.0)

    _rebin_spectrum(egrid->ener, xillverInputSpec, egrid->nbins, ener, spec_in, n_ener);

    // divide by the primary normalization factor to get the scaling of the xillver reflection spectrum
    double normFactorXill =
            calc_xillver_norm_from_primary_spectrum(xillverInputSpec, egrid->ener, egrid->nbins, status);

    delete[] xillverInputSpec;

    return normFactorXill;
}


/**
 * @brief Xillver reflection of the returning radiation of a zone in the disk frame, weighted over the inclinations of the
 * table with the weights @p dist, and normalized to the returning flux: N_xill(specRet) * sum_i dist_i I_i * |boost|.
 *
 * The weights decide whether the result is an intensity (GR emission-angle distribution, sum = 1) or a flux (slab flux
 * weights), see calc_zone_xillver_reflected_intensity() and calc_zone_xillver_angle_averaged_reflected_flux().
 */
static void calc_zone_xillver_reflected_weighted(const xillParam &xill_param,
                                                 double *dist,
                                                 const returnSpec2D *returnSpec,
                                                 VecD &xill_spec_returnrad,
                                                 int izone,
                                                 int *status) {
    get_angle_weighted_xill_spec(xill_spec_returnrad,
                                 returnSpec->ener,
                                 returnSpec->n_ener,
                                 &xill_param,
                                 dist,
                                 status);

    double xillverReflectionNormFactor
            = get_xillver_norm_factor_from_primary_spectrum(returnSpec->specRet[izone],
                                                            returnSpec->ener,
                                                            returnSpec->n_ener,
                                                            status);
    CHECK_STATUS_VOID(*status);

    for (int jj = 0; jj < returnSpec->n_ener; jj++) {
        xill_spec_returnrad[jj] *= fabs(xill_param.boost);
        xill_spec_returnrad[jj] *= xillverReflectionNormFactor;
    }
}

/**
 * Computes the reflected intensity of the returning radiation for a specific zone of a multi-zone accretion disk model.
 *
 * This function calculates the zone-specific reflected intensity in the disk frame, which is handed to the
 * relativistic convolution: the xillver table intensities, weighted with the GR emission-angle distribution of the
 * zone (rel_profile->rel_cosne->dist[izone], sum = 1, i.e., the angle-mixed intensity, which depends on the observer
 * inclination), and normalized by a normalization factor derived from the returning spectrum.
 * The final spectrum is multiplied with the absolute value of the boost parameter. The direct disk emission is not
 * included, as it is calculated for the whole disk at once (see calc_disk_emission_observed).
 *
 * Assumptions:
 * - A single blackbody temperature (kTbb) is assumed for the reflection of all zones.
 *
 * @see @ref returnrad_normalizations for a derivation of the normalization
 *      factors applied below and the rationale for splitting them this way, and
 *      @ref xillver_angular_dependence for intensity vs. flux.
 *
 * Parameters:
 * @param xill_param          Pointer to the structure containing xillver model parameters (e.g., spectral parameters and disk settings).
 * @param rel_profile         Pointer to the multi-zone diskline profile, including relativistic correction data for all zones.
 * @param returnSpec          Pointer to the 2D return spectrum structure containing zone-specific primary and reflected flux data.
 * @param xill_intensity      Array where the calculated reflected intensity for the specified zone is stored (length `returnSpec->n_ener`).
 * @param izone               Index of the zone for which the return flux is calculated.
 * @param status              Pointer to the status flag, used to report success or failure.
 */
void calc_zone_xillver_reflected_intensity(const xillParam &xill_param,
                                           const relline_spec_multizone *rel_profile,
                                           const returnSpec2D *returnSpec,
                                           VecD &xill_intensity,
                                           int izone,
                                           int *status) {
    calc_zone_xillver_reflected_weighted(xill_param,
                                         rel_profile->rel_cosne->dist[izone],
                                         returnSpec,
                                         xill_intensity,
                                         izone,
                                         status);
}

/**
 * @brief Angle-averaged reflected flux of the returning radiation of a zone in the disk frame (diagnostic output).
 *
 * As calc_zone_xillver_reflected_intensity(), but the xillver table intensities are weighted with the semi-infinite
 * slab flux weights 0.5*cos(incl)/n_incl (get_semi_infinite_slab_flux_weights()) instead of the GR emission-angle
 * distribution. The result is the flux the slab emits, independent of the observer, and can be compared directly with
 * the returning flux specRet (energy conservation) and with the mirrored flux. It is not used for the model spectrum.
 */
void calc_zone_xillver_angle_averaged_reflected_flux(const xillParam &xill_param,
                                                     const returnSpec2D *returnSpec,
                                                     VecD &xill_flux,
                                                     int izone,
                                                     int *status) {
    xillSpec *xill_spec = get_xillver_spectra(&xill_param, status);
    CHECK_STATUS_VOID(*status);
    VecD slab_weights = get_semi_infinite_slab_flux_weights(xill_spec);
    free_xill_spec(xill_spec);

    calc_zone_xillver_reflected_weighted(xill_param, slab_weights.data(), returnSpec, xill_flux, izone, status);
}

/**
   * @brief Compute the per-zone intensity of the returning radiation ("mirrored" instead of reflected).
   *
   * The returning flux F_ret (`specRet`) already contains the angles of incidence (each returning ray is weighted
   * with the cosine of its angle of incidence on the disk). "Mirrored" means that the emitted spectrum has the same
   * shape and flux as the incident one, assuming a reprocessing in the slab such that it is re-emitted with an
   * isotropic intensity (i.e., with the same cos(theta) dependence of the flux on the emission angle as in the
   * xillver definition). As for the xillver reflection (calc_zone_xillver_reflected_intensity), the spectrum
   * handed to the relativistic convolution must be an intensity in xillver table units, as the convolution converts
   * the intensity into the observed flux. The flux of an isotropic intensity I is 0.25*I
   * (SEMI_INFINITE_SLAB_FLUX_ISOTROPIC), so the mirrored intensity is I = F_ret / 0.25.
   *
   * The result is multiplied with the absolute value of the boost parameter. The direct disk emission is not
   * included, as it is calculated for the whole disk at once (see calc_disk_emission_observed).
   *
   * @see @ref xillver_angular_dependence for the conversion between intensity and flux.
   *
   * @param xill_param   Xillver parameters; only @c boost is read here.
   * @param returnSpec   Precomputed return-radiation spectra
   * @param mirrored_intensity  Output (intensity of the mirrored returning radiation) for zone ii, sized in returnSpec.
   * @param ii           Radial zone index.
   */
void calc_zone_mirrored_intensity(const xillParam &xill_param,
                                  const returnSpec2D *returnSpec,
                                  VecD &mirrored_intensity,
                                  int ii) {
    for (int jj = 0; jj < returnSpec->n_ener; jj++) {
        mirrored_intensity[jj] = returnSpec->specRet[ii][jj] / SEMI_INFINITE_SLAB_FLUX_ISOTROPIC * fabs(
                                     xill_param.boost);
    }
}


static void set_low_values_to_zero(double* spec, int n)
{
    double maxVal = 0.0;
    constexpr double lowestModelValue = 1e-8;

    for (int ii = 0; ii < n; ii++)
    {
        maxVal = fmax(maxVal, spec[ii]);
    }

    for (int ii = 0; ii < n; ii++)
    {
        if (spec[ii] < maxVal * lowestModelValue)
        {
            spec[ii] = 0.0;
        }
    }
}

static void set_values_outside_to_zero(double* spec, const double* ener, int n)
{
    for (int ii = 0; ii < n; ii++)
    {
        if (ener[ii] < EMIN_XILLVER || ener[ii + 1] > EMAX_XILLVER) {
            spec[ii] = 0.0;
        }
    }
}


void check_and_reset_value(double* ptr_val, double val_ref, const char* param_name)
{
    constexpr double prec = 1e-6;
    if (fabs(*ptr_val-val_ref) > prec)
    {
        *ptr_val = val_ref;
        printf(" re-setting parameter %s to %.3e (only this value allowed)\n", param_name, val_ref);
    }
}

void check_parameter_bounds_relxillbb(xillParam* xill_param,relParam* rel_param)
{
    if (rel_param->model_type == MOD_TYPE_RELXILLBB) {
        check_and_reset_value(&rel_param->a, 0.998, "a");
        check_and_reset_value(&rel_param->rin, kerr_rms(0.998), "Rin");
        check_and_reset_value(&rel_param->rout, 1000, "Rout");
        check_and_reset_value(&xill_param->lxi, 2, "logxi");
        check_and_reset_value(&xill_param->dens, 18, "logN");
    }
}

IonGradient calc_ion_gradient_from_returnrad(
    const returnSpec2D *returnSpec,
    const ModelDefinition &params,
    int *status) {
    const auto radialGrid = get_radial_grid_from_returntab(returnSpec);
    IonGradient ion_gradient{radialGrid, ION_GRAD_TYPE_ALPHA, 0};
    const auto primary_source_params = PrimarySourceParameters{params};
    auto emis_profile =
            calculate_emis_profile_from_spectra(returnSpec->specRet,
                                                returnSpec->rlo,
                                                returnSpec->rhi,
                                                returnSpec->nrad,
                                                returnSpec->ener,
                                                returnSpec->n_ener,
                                                status);
    ion_gradient.calculate_gradient(*emis_profile, primary_source_params);
    free_emisProfile(emis_profile);

    if (shouldOutfilesBeWritten()) {
        ion_gradient.write_to_file((fbase_rrad_bbody_dat + "emis_profile.dat").c_str());
    }

    return ion_gradient;
}


// boost=0: only the direct (primary) emission, without reflection
static bool is_boost_primary_only(double boost) {
    return fabs(boost) < 1e-8;
}

/**
 * @brief Diagnostic output of the reflected spectra of the zones (scaled with |boost|), in the rest frame
 * (rframe-reflect.fits) and observed, i.e., convolved with the line profile of each zone (obs-reflect.fits). If the
 * returning radiation is mirrored instead of reflected (mirrored_reflection_instead_of_xillver), the files are called
 * rframe-mirror-reflect.fits and obs-mirror-reflect.fits.
 *
 * All files are fluxes: in the rest frame the flux emitted by the zone (angle-averaged, independent of the observer),
 * observed the flux seen by the observer. The rest-frame spectra are therefore not the intensities that are convolved
 * (see @ref xillver_angular_dependence).
 */
static void write_bbret_diagnose_outfiles(const VecD_2D &specZonesConvOutput,
                                          const VecD_2D &specZonesReflectedFlux,
                                          const double *ener,
                                          int n_ener,
                                          const returnSpec2D *returnSpec,
                                          int *status) {
    const std::string name = mirrored_reflection_instead_of_xillver() ? "mirror-reflect" : "reflect";

    std::string fname = fbase_rrad_bbody + "obs-" + name + ".fits";
    std::cout << " writing BBret diagnose outfile: " << fname.c_str() << std::endl;
    fits_rr_write_2Dspec(fname.c_str(),
                         specZonesConvOutput,
                         ener,
                         n_ener,
                         returnSpec->rlo,
                         returnSpec->rhi,
                         returnSpec->nrad,
                         nullptr,
                         status);

    fname = fbase_rrad_bbody + "rframe-" + name + ".fits";
    fits_rr_write_2Dspec(fname.c_str(),
                         specZonesReflectedFlux,
                         ener,
                         n_ener,
                         returnSpec->rlo,
                         returnSpec->rhi,
                         returnSpec->nrad,
                         nullptr,
                         status);
}

/**
 * @brief Diagnostic output of the observed primary (direct) emission of each zone (obs-bbodyEmit.fits).
 *
 * In the model, the zones only contain the reflected emission, as the direct emission is calculated for the whole disk
 * at once (calc_disk_emission_observed). For the diagnostic output, the observed direct emission of each zone is its
 * primary (specPri, from which the returning radiation is calculated and which is written to rframe-bbodyEmit.fits)
 * convolved with the line profile of the zone.
 *
 * As for the reflected spectra, the convolution takes an intensity: specPri is the emitted flux of an isotropically
 * emitting blackbody surface, so its intensity specPri / 0.25 (SEMI_INFINITE_SLAB_FLUX_ISOTROPIC) is convolved, see
 * @ref xillver_angular_dependence. (The direct emission of the model, calc_disk_emission_observed, does not yet apply
 * this conversion and is therefore a factor 4 lower than the sum of this file over all zones.)
 */
static void write_primary_observed_diagnose_outfile(double *ener,
                                                    int n_ener,
                                                    const returnSpec2D *returnSpec,
                                                    const relline_spec_multizone *rel_profile,
                                                    specCache *spec_cache,
                                                    int *status) {
    VecD_2D spec_zones(returnSpec->nrad, VecD(n_ener));
    VecD intensity_pri(n_ener);
    for (int izone = 0; izone < returnSpec->nrad; izone++) {
        for (int jj = 0; jj < n_ener; jj++) {
            intensity_pri[jj] = returnSpec->specPri[izone][jj] / SEMI_INFINITE_SLAB_FLUX_ISOTROPIC;
        }
        convolveSpectrumFFTNormalized(ener,
                                      intensity_pri.data(),
                                      rel_profile->flux[izone],
                                      spec_zones[izone].data(),
                                      n_ener,
                                      1,
                                      1,
                                      izone,
                                      spec_cache,
                                      status);
    }

    const std::string fname = fbase_rrad_bbody + "obs-bbodyEmit.fits";
    std::cout << " writing BBret diagnose outfile: " << fname.c_str() << std::endl;
    fits_rr_write_2Dspec(fname.c_str(),
                         spec_zones,
                         ener,
                         n_ener,
                         returnSpec->rlo,
                         returnSpec->rhi,
                         returnSpec->nrad,
                         nullptr,
                         status);
}

static void set_xill_param_from_iongrad(xillParam &xill_param_zone, const IonGradient &ion_gradient, int izone) {
    xill_param_zone.dens = ion_gradient.dens[izone];
    xill_param_zone.lxi = ion_gradient.lxi[izone];
    xill_param_zone.radzone = static_cast<double>(izone);
}

/**
 * @brief Computes the relativistic reflection spectrum using a blackbody as the primary radiation source.
 *
 * This function generates the relativistic convolution of the incident spectrum using the parameters of the
 * blackbody source and the relativistic disk. It combines contributions from multiple zones to produce the
 * full output spectrum.
 *
 * @param spectrum The input/output XspecSpectrum containing the energy grid and flux array.
 * @param params Model definition holding all model parameters (rel + xill are extracted internally).
 * @param status Pointer to the status flag. Will be checked and updated during execution for error handling.
 *
 * The function performs the following major operations:
 * 1. Initialize and assert various parameters required for the computation.
 * 2. Generate a standard energy grid used for convolution.
 * 3. Load the XILLVER table data based on the model and primary source type.
 * 4. Calculate the ionization gradient from the return radiation emissivity profile.
 * 5. Compute the observed direct disk emission (boost >= 0) for the whole disk with a single convolution
 *    (calc_disk_emission_observed).
 * 6. Compute the reflected spectrum for each zone in the relativistic disk.
 * 7. Convolve the resulting spectrum using relativistic effects for each reflection zone.
 * 8. Rebin the convolved spectrum onto the input energy grid and combine contributions from all zones.
 * 9. Apply corrections to clean and format the output spectrum.
 *
 * Memory allocation and deallocation are handled for intermediate spectra within the function, ensuring proper cleanup
 * of resources to prevent memory leaks.
 *
 * Notes:
 * - The function assumes that the input and output spectra have been properly initialized by the caller.
 * - The input parameter xill_param->kTbb may be modified temporarily during computations but is restored
 *   to its original value upon completion.
 * - Diagnostic FITS files are written if diagnostic flags are set via shouldOutfilesBeWritten().
 */
void relxillbb_kernel(const XspecSpectrum &spectrum,
                      const ModelDefinition &params,
                      int* status)
{
    double *ener_inp = spectrum.energy;
    double *specOutput = spectrum.flux;
    int n_ener_inp = spectrum.num_flux_bins();

    auto rel_param = get_rel_params(params);
    auto xill_param = get_xill_params(params);

    check_parameter_bounds_relxillbb(xill_param, rel_param);

    CHECK_STATUS_VOID(*status);
    assert(xill_param->model_type == MOD_TYPE_RELXILLBB_NS || xill_param->model_type == MOD_TYPE_RELXILLBB);

    // get a standard grid for the convolution (is rebinned later to the input grid)
    auto ener_grid = get_relxill_conv_energy_grid();
    int n_ener = ener_grid->nbins;
    double* ener = ener_grid->ener;

    xillTable *xill_tab = nullptr;
    get_init_xillver_table(&xill_tab, xill_param->model_type, xill_param->prim_type, status);

    returnSpec2D* returnSpec = spec_returnrad_blackbody(ener, n_ener, xill_param->kTbb,
                                                        rel_param->rin, rel_param->rout, rel_param->a,
                                                        status);
    CHECK_RELXILL_ERROR("calculating return rad spec, failed", status);
    auto ion_gradient = calc_ion_gradient_from_returnrad(returnSpec, params, status);


    auto radialGrid = get_radial_grid_from_returntab(returnSpec);
    RelSysPar* systemParameters = get_system_parameters(rel_param, status);
    relline_spec_multizone
        * relProfile = relbase_profile(ener,
                                       n_ener,
                                       rel_param,
                                       systemParameters,
                                       xill_tab,
                                       radialGrid,
                                       status);


    // ========== //
    VecD specSingleZone(n_ener_inp);
    VecD_2D specZonesConvOutput(returnSpec->nrad);
    VecD_2D specZonesReflected(returnSpec->nrad); // intensity, handed to the convolution
    VecD_2D specZonesReflectedFlux(returnSpec->nrad); // emitted flux in the rest frame (only for the diagnostic output)
    // ========== //

    specCache *spec_cache = init_global_specCache(status);
    CHECK_STATUS_VOID(*status);
    setArrayToZero(specOutput, n_ener_inp);

    // direct disk emission (boost >= 0): a single convolution for the whole disk, with the temperature of each radial
    // point (the zones below only contain the reflected emission)
    if (xill_param->boost >= 0) {
        VecD specDirectObs = calc_disk_emission_observed(ener,
                                                         n_ener,
                                                         xill_param->kTbb,
                                                         rel_param,
                                                         systemParameters,
                                                         spec_cache,
                                                         status);
        _rebin_spectrum(ener_inp, specOutput, n_ener_inp, ener, specDirectObs.data(), n_ener);
    }

    for (int izone = 0; izone < relProfile->n_zones; izone++)
    {
        assert(returnSpec->n_ener == n_ener);

        // copy parameters for bug safety to not make changes to the global xill_param
        xillParam xill_param_zone = *xill_param;  // flat copy works, as no arrays contained

        specZonesReflected[izone] = VecD(n_ener, 0.0);
        specZonesReflectedFlux[izone] = VecD(n_ener, 0.0);
        specZonesConvOutput[izone] = VecD(n_ener, 0.0);

        // the zones only contain the reflected emission, i.e., nothing for boost=0 (the normalization of the FFT
        // convolution is not defined for a spectrum without flux)
        if (is_boost_primary_only(xill_param->boost)) {
            continue;
        }

        if (mirrored_reflection_instead_of_xillver()) {
            // TODO: remove this once everything is working; only for testing
            calc_zone_mirrored_intensity(xill_param_zone, returnSpec, specZonesReflected[izone], izone);
            if (shouldOutfilesBeWritten()) {
                for (int jj = 0; jj < n_ener; jj++) {
                    specZonesReflectedFlux[izone][jj] =
                            specZonesReflected[izone][jj] * SEMI_INFINITE_SLAB_FLUX_ISOTROPIC;
                }
            }
        } else {
            // very dangerous; need to check that it is within the table boundaries as well
            // (will only be used for the model with the xillverRTR table (MOD_TYPE_RELXILLBB))
            xill_param_zone.radzone = static_cast<double>(izone);
            if (use_iongrad_relxillbb()) {
                set_xill_param_from_iongrad(xill_param_zone, ion_gradient, izone);
            }
            calc_zone_xillver_reflected_intensity(xill_param_zone,
                                                  relProfile,
                                                  returnSpec,
                                                  specZonesReflected[izone],
                                                  izone,
                                                  status);
            if (shouldOutfilesBeWritten()) {
                calc_zone_xillver_angle_averaged_reflected_flux(xill_param_zone,
                                                                returnSpec,
                                                                specZonesReflectedFlux[izone],
                                                                izone,
                                                                status);
            }
        }


        convolveSpectrumFFTNormalized(ener,
                                      specZonesReflected[izone].data(),
                                      relProfile->flux[izone],
                                      specZonesConvOutput[izone].data(),
                                      n_ener,
                                      1,
                                      1,
                                      izone,
                                      spec_cache,
                                      status);


        _rebin_spectrum(ener_inp, specSingleZone.data(), n_ener_inp, ener, specZonesConvOutput[izone].data(), n_ener);

        // add spectra from the current zone to the output spectrum
        for (int jj = 0; jj < n_ener_inp; jj++) {
            specOutput[jj] += specSingleZone[jj];
        }
    }

    // clean spectrum
    set_low_values_to_zero(specOutput, n_ener_inp);
    set_values_outside_to_zero(specOutput, ener_inp, n_ener_inp);

    if (shouldOutfilesBeWritten()) {
        write_primary_observed_diagnose_outfile(ener, n_ener, returnSpec, relProfile, spec_cache, status);
        if (!is_boost_primary_only(xill_param->boost)) {
            write_bbret_diagnose_outfiles(specZonesConvOutput,
                                          specZonesReflectedFlux,
                                          ener,
                                          n_ener,
                                          returnSpec,
                                          status);
        }
    }

    delete_returnSpec2D(&returnSpec);
    delete_rel_params(rel_param);
    delete xill_param;
}

