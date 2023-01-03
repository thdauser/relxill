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

extern "C" {
#include "common.h"
}

#ifndef RELXILL__RELPHYSICS_H_
#define RELXILL__RELPHYSICS_H_

#define TPROFILE_ALPHA 1
#define TPROFILE_DISKBB 2

double ut_disk(double r, double a);
double calc_proper_area_ring(double rlo, double rhi, double spin);

void disk_Tprofile_alpha(const double *rlo, const double *rhi, double *temp, int n, double Rin, double Tin);
void disk_Tprofile_diskbb(const double *rlo, const double *rhi, double *temp, int n, double Tin);
double *get_tprofile(double *rlo, double *rhi, int nrad, double Rin, double Tin, int type, int *status);

void bbody_spec(const double *ener, int n, double *spec, double temperature, double gfac);

double density_ss73_zone_a(double radius, double rms);

// energy shift from the primary source to the observer
double calc_g_inf(double height, double a);

/* calculate the radius of marginal stability */
double kerr_rms(double a);

/* get the rplus value (size if the black hole event horizon) */
double kerr_rplus(double a);

/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet);

/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del);

double relat_abberation(double del, double beta);

double doppler_factor_source_obs(const relParam *rel_param);
double energy_shift_source_obs(const relParam *rel_param);
double energy_shift_source_disk(const relParam *rel_param, double radius_disk, double del_emit);

double calc_lp_emissivity_newton(double h, double r);

#endif //RELXILL__RELPHYSICS_H_
