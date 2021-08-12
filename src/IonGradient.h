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
#ifndef IONGRADIENT_H_
#define IONGRADIENT_H_

extern "C" {
#include "relutility.h"
#include "relphysics.h"
}

typedef struct {
  double *lxi;
  double *fx;
  double *r;
  double *del_emit;
  double *dens;
  int nbins;
} ion_grad;


ion_grad *calc_ion_gradient(relParam *rel_param,
                            double xlxi0,
                            double xindex,
                            int type,
                            double *rgrid,
                            int n,
                            int *status);

void free_ion_grad(ion_grad *ion);

#endif
