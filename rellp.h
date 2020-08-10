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
#ifndef RELLP_H_
#define RELLP_H_

#include "relbase.h"
#include "relutility.h"
#include "common.h"

#define NHBINS_VERTICALLY_EXTENDED_SOURCE 50



typedef struct{

  double* heightArr;  //nh+1 bins
  double* heightMean;
  double* beta;
  int nh;

} extPrimSource;


// calculate the angles of emission from the primary source to git Rin and Rout
// void get_ad_del_lim(relParam* param, relSysPar* sysPar, int* status);

emisProfile* calc_emis_profile(double* r, int nr, relParam* param, int* status);

void get_emis_jet(emisProfile*, relParam* param, int* status);

int modelLampPostPointsource(relParam* param);

extPrimSource* getExtendedJetGeom(const relParam *param, int* status);

////////////

void free_cached_lpTable(void);

lpReflFrac* new_lpReflFrac(int* status);
void free_lpReflFrac(lpReflFrac** str);

emisProfile* new_emisProfile(double* re, int nr, int* status);
void free_emisProfile(emisProfile* emis_profile);

extPrimSource* new_extendedPrimarySource(int nh, int* status);
void free_extendedPrimarySource (extPrimSource* source);



#endif /* RELLP_H_ */
