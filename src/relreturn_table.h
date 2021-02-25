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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELRETURN_TABLE_H_
#define RELRETURN_TABLE_H_

#include "relutility.h"
#include "relphysics.h"

#define RMAX_RELRET 1000

typedef struct {

  double a;  // store the spin here as well ?!

  double *rlo;      // grid for r_i and r_e
  double *rhi;
  int nrad;

  /*  2d array of fractions, first dimension is r_i, second is r_e */
  double **frac_e;  // not necessary, only use for debugging
  double **frac_i;  //

  /* extending to the 3rd dimension to g */
  double ***frac_g;
  int ng;

  double **gmin;
  double **gmax;

  /* return fraction, depending on r_e */
  double *f_ret;
  double *f_inf;
  double *f_bh;

} tabulatedReturnFractions;

typedef struct {

  double *spin;
  int nspin;

  tabulatedReturnFractions **retFrac;

} returnTable;


typedef struct {
  double a;  // store the spin here as well ?!

  double *rlo; // for debugging only
  double *rhi; // for debugging only
  double* rad; // grid for r_i and r_e
  int nrad;

  double *proper_area_ring;  // proper area of this ring

  /* indicies of the radial bins that are used here from tabData, i.e., that are
   * now the new rlo, rhi and therefore rlo = tabData->rlo[irad] (except for Rin and Rout bin)*/
  int *irad;

  /*  2d array of fractions, first dimension is r_incident, second is r_emitted */
  double **frac_i;  //

  double *f_ret; // fraction of returning photons  [dimension r_incident]
  double *f_inf; // fraction of photons reaching the observed [dimension r_incident]

  /* everything we do not interpolate we directly take from the table */
  tabulatedReturnFractions *tabData;

} returningFractions;


/* Routines */

returnTable *get_returnrad_table(int *status);

void free_2d(double ***vals, int n1);
void free_cached_returnTable(void);
void free_returningFractions(returningFractions **dat);

returningFractions *get_rrad_fractions(double spin, double rin, double rout, int *status);

#endif /* RELRETURN_TABLE_H_ */
