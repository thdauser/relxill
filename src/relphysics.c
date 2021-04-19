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

    Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "relphysics.h"

#include "common.h"
#include "relutility.h"


double calc_proper_area_ring(double rlo, double rhi, double a) {

  double rmean = 0.5 * (rlo + rhi);

  assert(a < 1);
  assert(a >= -1);

  double rho2 = rmean * rmean + a * a;  // this is rho^2
  double del = rmean * rmean - 2 * rmean + a * a;

  double area_gr = 2 * M_PI / sqrt(del);
  area_gr *= sqrt((rho2 * rho2 + 2 * a * a * rmean));

  area_gr *= (rhi - rlo);

  return area_gr;
}

// Black Body Spectrum (returns ph/s/cm²/keV)  (verified with Xspec bbody model)
// ToDO: use a properly integrated model (?)x
void bbody_spec(double *ener, int n, double *spec, double temperature, double gfac) {
  for (int ii = 0; ii < n; ii++) {
    double emean = 0.5 * (ener[ii] + ener[ii + 1]);
    spec[ii] = pow((emean / gfac), 2) / (exp(emean / gfac / temperature) - 1);
  }
}

// calculate the temperature for a given radius r and respective Rin following SS73
static double disk_temperature_alpha(double r, double Rin) {
  return pow((r / Rin), (-3. / 4)) * pow((1 - sqrt(Rin / r)), (1. / 4));
}

// temperature of a disk as used by diskbb (equal to alpha disk for large radii)
static double disk_temperature_diskbb(double r, double Rin, double Tin) {
  return pow(Tin * (r / Rin), (-3. / 4));  // if we use Rin, we don't get T=Tin at the inner edge
}

// temperature profile for the given radii; use Tmax at Rmax following Poutanan+2007
void disk_Tprofile_alpha(const double *rlo, const double *rhi, double *temp, int n, double Rin, double Tin) {

  // use Tin as Tmax
  double Rmax = pow((1.5), (4. / 5)) * Rin;  // Poutanen+2007
  double K = Tin / disk_temperature_alpha(Rmax, Rin);

  double RminGrid = 0.5 * (rlo[0] + rhi[0]);
  assert(RminGrid >= Rin);

  for (int ii = 0; ii < n; ii++) {
    temp[ii] = K * disk_temperature_alpha(0.5 * (rlo[ii] + rhi[ii]), Rin);
  }
}

// temperature profile for the given radii; use Tmax at Rmax following Poutanan+2007
void disk_Tprofile_diskbb(const double *rlo, const double *rhi, double *temp, int n, double Tin) {

  double Rin = 0.5 * (rlo[0] + rhi[0]); // needs to be rmean, such that we really get Tin at the inner zone
  for (int ii = 0; ii < n; ii++) {
    temp[ii] = disk_temperature_diskbb(0.5 * (rlo[ii] + rhi[ii]), Rin, Tin);
  }
}

double *get_tprofile(double *rlo, double *rhi, const int nrad, double Rin, double Tin, int type, int *status) {

  CHECK_STATUS_RET(*status, NULL);

  double *tprofile = (double *) malloc(sizeof(double) * nrad);
  CHECK_MALLOC_RET_STATUS(tprofile, status, NULL)

  if (type == TPROFILE_ALPHA) {
    disk_Tprofile_alpha(rlo, rhi, tprofile, nrad, Rin, Tin);
  } else if (type == TPROFILE_DISKBB) {
    disk_Tprofile_diskbb(rlo, rhi, tprofile, nrad, Tin);
  } else {
    RELXILL_ERROR(" failed getting the temperature profile ", status);
    printf("    reason: type=%i is unknown\n", type);
  }

  return tprofile;
}
