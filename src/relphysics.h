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

#ifndef RELXILL__RELPHYSICS_H_
#define RELXILL__RELPHYSICS_H_

#define TPROFILE_ALPHA 1
#define TPROFILE_DISKBB 2

double ut_disk(double r, double a);
double calc_proper_area_ring(double rlo, double rhi, double spin);

void disk_Tprofile_alpha(const double *rlo, const double *rhi, double *temp, int n, double Rin, double Tin);
void disk_Tprofile_diskbb(const double *rlo, const double *rhi, double *temp, int n, double Tin);
double *get_tprofile(double *rlo, double *rhi, const int nrad, double Rin, double Tin, int type, int *status);

void bbody_spec(double *ener, int n, double *spec, double temperature, double gfac);

double density_ss73_zone_a(double radius, double rms);

double calc_g_inf(double height, double a);

#endif //RELXILL__RELPHYSICS_H_
