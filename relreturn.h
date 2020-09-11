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

#ifndef RELXILL__RELRETURN_H_
#define RELXILL__RELRETURN_H_

#include "relreturn_table.h"
#include "relreturn_datastruct.h"

typedef struct {

  double *ener;  // n+1 elements
  int n;

  double *rlo;
  double *rhi;
  int nrad;

  returnFracIpol *rDat;

  double **specs; // n elements

} spec2D;


void fits_rr_write_2Dspec(char *fname,
                          double **spec_arr,
                          double *ener,
                          int nener,
                          double *rlo,
                          double *rhi,
                          int nrad,
                          returnFracIpol *dat,
                          int *status);

/**  return the diskbb spectrum for the radial grid given by the table for this spin value **/
void spec_diskbb(double *ener, double *spec, int n, double Tin, double spin, int *status);


// spec_prim can be NULL
returnSpec2D *spec_returnrad_blackbody(double *ener, double *spec, double *spec_prim, int nener,
                                       double Tin, double Rin, double Rout, double spin, int *status);


void relxill_bb_kernel(double *ener_inp,
                       double *spec_inp,
                       int n_ener_inp,
                       xillParam *xill_param,
                       relParam *rel_param,
                       int noXillverReflection,
                       int *status);

double *getRadialGridFromReturntab(returnSpec2D *spec, int* status);

double * getXillverPrimaryBBodyNormalized(double kTbb, double* spec_in, double* ener, int n_ener, int* status);

double calcXillverNormfacRetrad2BoodyAtHighenergy(double kTbb, double * spec_in, double* spec_bb, double* ener, int n_ener, int* status);

double *scaledXillverPrimaryBBodyHighener(double kTbb, double *spec_in, double *ener, int n_ener, int *status);

void getZoneReflectedReturnFluxDiskframe(xillParam *xill_param, rel_spec *rel_profile, const returnSpec2D *returnSpec,
                                         double *xill_flux_returnrad, int izone, int *status);

void getZoneIncidentReturnFlux(xillParam *xill_param, const returnSpec2D *returnSpec, double *returnFlux, int ii);
void getZoneDirectPrimaryFlux(xillParam *xill_param, const returnSpec2D *returnSpec, double *returnFlux, int ii);

double *getTemperatureProfileDiskZones(returnFracIpol *dat, double Rin, double Tin, int *status);

void init_par_relxill_bbret(relParam **rel_param,
                            xillParam **xill_param,
                            const double *inp_par,
                            int n_parameter,
                            int *status);
void tdrelxillbbret(const double *ener0,
                    int n_ener0,
                    double *photar,
                    const double *parameter,
                    int n_parameter,
                    int *status);

void interpolEmisProfile(emisProfile *emisReb, emisProfile *emis0, int *status);

#endif //RELXILL__RELRETURN_H_
