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

/**  return the diskbb spectrum for the radial grid given by the table for this spin value **/
void spec_diskbb(double *ener, double *spec, int n, double Tin, double spin, int *status);


// spec_prim can be NULL
returnSpec2D *spec_returnrad_blackbody(double *ener,
                                       int nener,
                                       double Tin,
                                       double Rin,
                                       double Rout,
                                       double spin,
                                       int *status);


void relxill_bb_kernel(const XspecSpectrum &spectrum,
                       const ModelDefinition &params,
                       int *status);

std::vector<double> getRadialGridFromReturntab(const returnSpec2D *spec);

VecD getXillverPrimaryBBody(double kTbb, double *spec_in, const double *ener, int n_ener, int *status);

double calcNormfacBBodyAtHighenergy(double kTbb, double *spec_bb, const double *ener, int n_ener, int *status);

VecD getXillverPimaryBBodyNormalizedAtHighener(double kTbb,
                                               double *spec_in,
                                               const double *ener,
                                               int n_ener,
                                               int *status);

void calcZoneDirectAndReflectedReturnFlux(xillParam *xill_param,
                                          relline_spec_multizone *rel_profile,
                                          const returnSpec2D *returnSpec,
                                          VecD &xill_flux_returnrad,
                                          int izone,
                                          int *status);

void calcZoneDirectAndReturnFlux(xillParam *xill_param, const returnSpec2D *returnSpec, VecD &returnFlux, int ii);

VecD getTemperatureProfileDiskZones(returningFractions *dat, double Rin, double Tin, int *status);

IonGradient calc_ion_gradient_from_returnrad(const returnSpec2D *returnSpec,
                                             const ModelDefinition &params,
                                             int *status);

#endif //RELXILL__RELRETURN_H_