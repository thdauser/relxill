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
#ifndef MODELS_H_
#define MODELS_H_

#include "relbase.h"
#include "relutility.h"
#include "xilltable.h"

/**** DEFINES **/
#define MOD_TYPE_RELLINE 1
#define NUM_PARAM_RELLINE 10

#define MOD_TYPE_RELLINELP 2
#define NUM_PARAM_RELLINELP 9

#define MOD_TYPE_RELCONV 11
#define NUM_PARAM_RELCONV 8

#define MOD_TYPE_RELCONVLP 12
#define NUM_PARAM_RELCONVLP 7

#define MOD_TYPE_XILLVER 0
#define NUM_PARAM_XILLVER 7

#define MOD_TYPE_XILLVER_NTHCOMP 100
#define NUM_PARAM_XILLVER_NTHCOMP 7

#define MOD_TYPE_RELXILL -1
#define NUM_PARAM_RELXILL 13

#define MOD_TYPE_RELXILLLP -2
#define NUM_PARAM_RELXILLLP 12

/** density models **/
#define MOD_TYPE_RELXILLDENS -10
#define NUM_PARAM_RELXILLDENS 13

#define MOD_TYPE_RELXILLLPDENS -11
#define NUM_PARAM_RELXILLLPDENS 12

#define MOD_TYPE_XILLVERDENS -100
#define NUM_PARAM_XILLVERDENS 7

/** ion grad models **/
#define MOD_TYPE_RELXILLLPION -21
#define NUM_PARAM_RELXILLLPION 15

// TODO: implement RELXILLION model ??

/** CO models **/
#define MOD_TYPE_RELXILLCO -200
#define NUM_PARAM_RELXILLCO 14

#define MOD_TYPE_XILLVERCO -210
#define NUM_PARAM_XILLVERCO 8

/** Neutron Star / BB models **/
#define MOD_TYPE_RELXILLNS -30
#define NUM_PARAM_RELXILLNS 13

#define MOD_TYPE_XILLVERNS -101
#define NUM_PARAM_XILLVERNS 7

#define MOD_TYPE_XILLVERDENS_NTHCOMP 1000
#define NUM_PARAM_XILLVERDENS_NTHCOMP 8

#define MOD_TYPE_RELXILLDENS_NTHCOMP -1001
#define NUM_PARAM_RELXILLDENS_NTHCOMP 14

#define MOD_TYPE_RELXILLLPDENS_NTHCOMP -1002

#define NUM_PARAM_RELXILLLPDENS_NTHCOMP 15

// unpublished models
#define MOD_TYPE_RELXILLBBRET -300
#define MOD_TYPE_RELXILLLPRET -310

/****  TYPE DEFINITIONS ****/


/**** FUNCTION DEFINITIONS ****/

void check_parameter_bounds(relParam *param, int *status);
double *shift_energ_spec_1keV(const double *ener, const int n_ener, double line_energ, double z, int *status);

relParam *init_par_relline(const double *inp_par, const int n_parameter, int *status);
relParam *init_par_relline_lp(const double *inp_par, const int n_parameter, int *status);
relParam *init_par_relconv(const double *inp_par, const int n_parameter, int *status);
xillParam *init_par_xillver(const double *inp_par, const int n_parameter, int *status);
xillParam *init_par_xillver_nthcomp(const double *inp_par, const int n_parameter, int *status);
xillParam *init_par_xillver_ns(const double *inp_par, const int n_parameter, int *status);
xillParam *init_par_xillver_co(const double *inp_par, const int n_parameter, int *status);
xillParam *init_par_xillver_dens(const double *inp_par, const int n_parameter, int *status);

void init_par_relxill(relParam **rel_param,
                      xillParam **xill_param,
                      const double *inp_par,
                      const int n_parameter,
                      int *status);

void init_par_relxilllp(relParam **rel_param,
                        xillParam **xill_param,
                        const double *inp_par,
                        const int n_parameter,
                        int *status);

/** basic xillver model function **/
void xillver_base(const double *ener0, const int n_ener0, double *photar, xillParam *param_struct, int *status);

void relline_base(double *ener1keV, double *photar, const int n_ener, relParam *param_struct, int *status);

/** internal MODEL FUNCTIONS **/
void tdrelline(const double *ener,
               const int n_ener,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status);
void tdrellinelp(const double *ener,
                 const int n_ener,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status);
void tdrelxill(const double *ener0,
               const int n_ener0,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status);
void tdrelxilllp(const double *ener0,
                 const int n_ener0,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status);
void tdrelxilllpion(const double *ener0,
                    const int n_ener0,
                    double *photar,
                    const double *parameter,
                    const int n_parameter,
                    int *status);
void tdxillver(const double *ener0,
               const int n_ener0,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status);
void tdrelconv(const double *ener,
               const int n_ener,
               double *photar,
               const double *parameter,
               const int n_parameter,
               int *status);
void tdrelconvlp(const double *ener,
                 const int n_ener,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status);
void tdrelxilldens(const double *ener0,
                   const int n_ener0,
                   double *photar,
                   const double *parameter,
                   const int n_parameter,
                   int *status);
void tdrelxilllpdens(const double *ener0,
                     const int n_ener0,
                     double *photar,
                     const double *parameter,
                     const int n_parameter,
                     int *status);
void tdxillverdens(const double *ener0,
                   const int n_ener0,
                   double *photar,
                   const double *parameter,
                   const int n_parameter,
                   int *status);
void tdrelxillns(const double *ener0,
                 const int n_ener0,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status);
void tdxillverns(const double *ener0,
                 const int n_ener0,
                 double *photar,
                 const double *parameter,
                 const int n_parameter,
                 int *status);

void tdrelxill_nthcomp(const double *ener0,
                       const int n_ener0,
                       double *photar,
                       const double *parameter,
                       const int n_parameter,
                       int *status);
void tdrelxilllp_nthcomp(const double *ener0,
                         const int n_ener0,
                         double *photar,
                         const double *parameter,
                         const int n_parameter,
                         int *status);
void tdrelxilllpion_nthcomp(const double *ener0,
                            const int n_ener0,
                            double *photar,
                            const double *parameter,
                            const int n_parameter,
                            int *status);
void tdxillver_nthcomp(const double *ener0,
                       const int n_ener0,
                       double *photar,
                       const double *parameter,
                       const int n_parameter,
                       int *status);

void tdxillverco(const double *ener0, const int n_ener0, double *photar, const double *parameter, const int n_parameter,
                 int *status);

void tdrelxillco(const double *ener0, const int n_ener0, double *photar, const double *parameter, const int n_parameter,
                 int *status);

// Dens & Nthcomp Model
xillParam *init_par_xillver_dens_nthcomp(const double *inp_par, const int n_parameter, int *status);
void init_par_relxilldens_nthcomp(relParam **rel_param,
                                  xillParam **xill_param,
                                  const double *inp_par,
                                  const int n_parameter,
                                  int *status);
void init_par_relxilllp_dens_nthcomp(relParam **rel_param,
                                     xillParam **xill_param,
                                     const double *inp_par,
                                     const int n_parameter,
                                     int *status);
void tdrelxilldens_nthcomp(const double *ener0,
                           const int n_ener0,
                           double *photar,
                           const double *parameter,
                           const int n_parameter,
                           int *status);
void tdrelxilllpdens_nthcomp(const double *ener0,
                             const int n_ener0,
                             double *photar,
                             const double *parameter,
                             const int n_parameter,
                             int *status);
void tdxillverdens_nthcomp(const double *ener0,
                           const int n_ener0,
                           double *photar,
                           const double *parameter,
                           const int n_parameter,
                           int *status);

/* get the version number text on the screen (if not already printed before */
void print_version_number(int *status);

/* get a new relbase parameter structure and initialize it */
relParam *new_relParam(int model_type, int emis_type, int *status);

/* free relbase parameter */
void free_relParam(relParam *);

xillParam *new_xillParam(int model_type, int prim_type, int *status);
void free_xillParam(xillParam *);




#endif /* MODELS_H_ */
