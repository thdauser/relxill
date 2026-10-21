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

#ifndef RELXILL__RELRETURN_H_
#define RELXILL__RELRETURN_H_

#include "Types.h"
#include "Relreturn_Table.h"
#include "Relreturn_Datastruct.h"
#include "Relbase.h"
#include "XspecSpectrum.h"
#include "ModelDefinition.h"
#include "IonGradient.h"


extern "C" {
#include "common.h"
}

#include <vector>

// maximal ratio r_hi/r_lo of the sub-zones on which the direct disk emission of each radial zone is calculated (see
// get_disk_emission_zones_binned): the deviation from a fine radial grid is 0.02-0.03% for the alpha-disk and the
// Novikov-Thorne profile (0.2% for the diskbb profile), see doc/returnrad_disk_emission.md
#define RRAD_BBODY_MAX_RATIO_SUBZONES 1.1

typedef struct {
  double *ener; // n+1 elements
  int n;

  double *rlo;
  double *rhi;
  int nrad;

  returningFractions *rDat;

  double **specs; // n elements
} spec2D;


void fits_rr_write_2Dspec(const char *fname,
                          double **spec_arr,
                          const double *ener,
                          int nener,
                          const VecD &rlo,
                          const VecD &rhi,
                          int nrad,
                          returningFractions *dat,
                          int *status);

/**  direct disk emission of the radial zones [rlo, rhi], calculated on a fine radial grid (sub-zones) within each
     zone, see the documentation in Relreturn_BlackBody.cpp **/
VecD_2D get_disk_emission_zones_binned(const double *ener,
                                       int nener,
                                       const VecD &rlo,
                                       const VecD &rhi,
                                       double Rin,
                                       double kT,
                                       double spin,
                                       int tprofile,
                                       int *status);

/**  returning radiation [cts/bin on ener_inp] of the primary spectra spec_pri of the zones (on the bbody energy grid,
     get_std_bbody_energy_grid), see the documentation in Relreturn_BlackBody.cpp **/
VecD_2D get_returnrad_specs(double *ener_inp,
                            int nener_inp,
                            const returningFractions *dat,
                            const VecD_2D &spec_pri,
                            int norder,
                            int *status);

/**  multicolor disk spectrum (no relativistic effects) for the radial zones [rlo, rhi] and the temperature profile
     tprofile (TPROFILE_ALPHA or TPROFILE_DISKBB), see the documentation in Relreturn_BlackBody.cpp **/
void spec_multicolor_disk(const double *ener,
                          double *spec,
                          int n,
                          const VecD &rlo,
                          const VecD &rhi,
                          double kT,
                          int tprofile,
                          int *status);


// spec_prim can be NULL
returnSpec2D *spec_returnrad_blackbody(double *ener,
                                       int nener,
                                       double Tin,
                                       double Rin,
                                       double Rout,
                                       double spin,
                                       int *status);


void relxillbb_kernel(const XspecSpectrum &spectrum,
                      const ModelDefinition &params,
                      int *status);

std::vector<double> get_radial_grid_from_returntab(const returnSpec2D *spec);

void calc_zone_xillver_reflected_flux(const xillParam &xill_param,
                                      const relline_spec_multizone *rel_profile,
                                      const returnSpec2D *returnSpec,
                                      VecD &xill_flux_returnrad,
                                      int izone,
                                      int *status);

void calc_zone_mirrored_flux(const xillParam &xill_param,
                             const returnSpec2D *returnSpec,
                             VecD &returnFlux,
                             int ii);

VecD get_temperature_profile_disk_zones(const returningFractions *dat, double Rin, double Tin, int *status);

IonGradient calc_ion_gradient_from_returnrad(const returnSpec2D *returnSpec,
                                             const ModelDefinition &params,
                                             int *status);

#endif //RELXILL__RELRETURN_H_