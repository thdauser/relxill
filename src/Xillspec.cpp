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

#include "Xillspec.h"
#include "Relphysics.h"
#include "Relxill_errorhandling.h"

extern "C" {
#include "xilltable_reading.h"
#include "common.h"
}

#include <vector>
#include <memory>
#include <cassert>


EnerGrid *global_xill_egrid_coarse = nullptr;
EnerGrid *global_energy_grid_relxill = nullptr;


/** @brief main routine for the xillver table: returns a spectrum for the given parameters
 *  @details
 *   - decides if the table needs to be initialized and/or more data loaded
 *   - automatically normalizes  the spectra to logN=1e15 and logXi=0
 *   - for relxill type models it returns spectra for all available inclinations
 * @param param
 * @param status
 * @return xillSpec
 */



EnerGrid *get_relxill_conv_energy_grid() {
  if (global_energy_grid_relxill == nullptr) {
    global_energy_grid_relxill = new EnerGrid;
    global_energy_grid_relxill->nbins = N_ENER_CONV;
    global_energy_grid_relxill->ener = new double[global_energy_grid_relxill->nbins + 1];
    get_log_grid(global_energy_grid_relxill->ener, global_energy_grid_relxill->nbins + 1,
                 EMIN_RELXILL_CONV, EMAX_RELXILL_CONV);
  }

  return global_energy_grid_relxill;
}


xillSpec* new_xill_spec(int n_incl, int n_ener, int* status)
{
  CHECK_STATUS_RET(*status, NULL);

  xillSpec* spec = (xillSpec*)malloc(sizeof(xillSpec));
  CHECK_MALLOC_RET_STATUS(spec, status, NULL)

  spec->n_ener = n_ener;
  spec->n_incl = n_incl;
  spec->ener = NULL;
  spec->incl = NULL;
  spec->flu = NULL;

  do
  {
    spec->ener = (double*)malloc(sizeof(double) * (n_ener + 1));
    CHECK_MALLOC_BREAK_STATUS(spec->ener, status)

    spec->incl = (double*)malloc(sizeof(double) * (n_incl));
    CHECK_MALLOC_BREAK_STATUS(spec->incl, status)

    spec->flu = (double**)malloc(sizeof(double*) * n_incl);
    CHECK_MALLOC_BREAK_STATUS(spec->flu, status)

    int ii;
    for (ii = 0; ii < n_incl; ii++)
    {
      spec->flu[ii] = (double*)malloc(sizeof(double) * n_ener);
      CHECK_MALLOC_BREAK_STATUS(spec->flu[ii], status)
    }
  }
  while (0);

  if (*status == EXIT_FAILURE || spec->ener == NULL || spec->incl == NULL || spec->flu == NULL)
  {
    printf(" error: allocating memeory for the xillver spectrum failed!\n");
    free_xill_spec(spec);
    return NULL;
  }

  assert(spec != NULL);

  return spec;
}

void free_xill_spec(xillSpec* spec)
{
  if (spec != NULL)
  {
    free(spec->ener);
    free(spec->incl);
    if (spec->flu != NULL)
    {
      int ii;
      for (ii = 0; ii < spec->n_incl; ii++)
      {
        free(spec->flu[ii]);
      }
      free(spec->flu);
    }
    free(spec);
  }
}


std::vector<int> get_xilltab_indices_for_paramvals(const xillTableParam* param, const xillTable* tab, int* status)
{
  CHECK_STATUS_RET(*status, std::vector<int>{});

  // Use unique_ptr for automatic memory management
  std::unique_ptr<float[], decltype(&free)> inp_global_paramvals(
    get_xilltab_global_paramvals(param, status),
    free
  );

  if (*status != EXIT_SUCCESS || !inp_global_paramvals)
  {
    return std::vector<int>{};
  }

  // Use vector instead of manual memory management
  std::vector<int> ind_xilltab(tab->num_param);

  for (int ii = 0; ii < tab->num_param; ii++)
  {
    assert(tab->num_param_vals[ii] > 1); // we need at least two values per parameter (for the interpolation later)

    float pval = inp_global_paramvals[tab->param_index[ii]];
    ind_xilltab[ii] = binary_search_float(tab->param_vals[ii], tab->num_param_vals[ii], pval);

    // routine should be safe, meaning index always within the array bounds
    assert(ind_xilltab[ii] >= 0 && ind_xilltab[ii] < tab->num_param_vals[ii]);

    // make sure all parameters are by default within the defined limits here!!
    if (tab->param_index[ii] == PARAM_RID && pval < 0)
    {
      RELXILL_ERROR("radial zone index not set (needs to be given for the xillver-RTR table)", status);
      return std::vector<int>{};
    }

    // for the other parameters, we reset it
    if (ind_xilltab[ii] < 0)
    {
      ind_xilltab[ii] = 0;
    }
    else if (ind_xilltab[ii] > tab->num_param_vals[ii] - 2)
    {
      ind_xilltab[ii] = tab->num_param_vals[ii] - 2;
    }
  }

  return ind_xilltab;
}


// Default values from the param-structure, as those are set to the value the table is calculated
// with, even if the model does not have a lxi or density parameter
static double getDefaultDensity(const xillTableParam* param)
{
  return param->dens;
}

static double getDefaultLogxi(const xillTableParam* param)
{
  return param->lxi;
}


void check_xilltab_cache(const char* fname, const xillTableParam* param, const xillTable* tab,
                         std::vector<int> ind, int* status)
{
  CHECK_STATUS_VOID(*status);

  fitsfile* fptr = nullptr;

  // Determine the actual number of parameters we need to iterate over
  int num_param = tab->num_param;
  assert(tab->num_param == ind.size()); // check that we have the correct number of indizes
  int max_dims = 6; // Maximum supported dimensions

  // Set up loop bounds for each dimension
  int loop_start[max_dims];
  int loop_end[max_dims];

  // Set up bounds for the all dimensions (not we count backwards, meaning first entries are zero for dim<max_dims)
  for (int iloop = 0; iloop < max_dims; iloop++)
  {
    int ipar = num_param - max_dims + iloop; // parameter index from the input array
    if (iloop == max_dims - 1) // the last entry is always the inclination
    {
      assert(ipar == num_param - 1);
      // always load **all** incl bins as for relxill we will certainly need it
      loop_start[iloop] = 0;
      loop_end[iloop] = tab->n_incl;
    }
    else if (ipar >= 0 && ipar < num_param - 1) // all remaining existing dimensions
    {
      loop_start[iloop] = ind[ipar];
      loop_end[iloop] = ind[ipar] + 1;
    }
    else
    {
      // dimension does not exist
      loop_start[iloop] = 0;
      loop_end[iloop] = 0;
    }
  }

  double defDensity = getDefaultDensity(param);
  double defLogxi = getDefaultLogxi(param);

  int ii, jj, kk, ll, mm, nn;
  for (nn = loop_start[0]; nn <= loop_end[0]; nn++)
  {
    for (ii = loop_start[1]; ii <= loop_end[1]; ii++)
    {
      for (jj = loop_start[2]; jj <= loop_end[2]; jj++)
      {
        for (kk = loop_start[3]; kk <= loop_end[3]; kk++)
        {
          for (ll = loop_start[4]; ll <= loop_end[4]; ll++)
          {
            for (mm = loop_start[5]; mm < loop_end[5]; mm++)
            {
              if (get_xillspec_table(tab, nn, ii, jj, kk, ll, mm) == nullptr)
              {
                xilltable_fits_load_single_spec(fname,
                                                &fptr,
                                                tab,
                                                defDensity,
                                                defLogxi,
                                                nn,
                                                ii,
                                                jj,
                                                kk,
                                                ll,
                                                mm,
                                                status);
                CHECK_STATUS_VOID(*status);
              }
            }
          }
        }
      }
    }
  }

  CHECK_RELXILL_DEFAULT_ERROR(status);

  if (fptr != nullptr)
  {
    if (fits_close_file(fptr, status))
    {
      RELXILL_ERROR(" *** error closing FITS file", status);
    }
  }

  assert(tab->data_storage != nullptr);
}


static void interp_5d_tab_incl(const xillTable* tab, double* flu, int n_ener,
                               double f1, double f2, double f3, double f4,
                               int i0, int i1, int i2, int i3, int i4, int i5)
{
  float* dat_0000 = get_xillspec_table(tab, i0, i1, i2, i3, i4, i5);
  float* dat_1000 = get_xillspec_table(tab, i0, i1 + 1, i2, i3, i4, i5);
  float* dat_0100 = get_xillspec_table(tab, i0, i1, i2 + 1, i3, i4, i5);
  float* dat_0010 = get_xillspec_table(tab, i0, i1, i2, i3 + 1, i4, i5);

  float* dat_1100 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3, i4, i5);
  float* dat_1010 = get_xillspec_table(tab, i0, i1 + 1, i2, i3 + 1, i4, i5);
  float* dat_0110 = get_xillspec_table(tab, i0, i1, i2 + 1, i3 + 1, i4, i5);
  float* dat_1110 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4, i5);

  float* dat_0001 = get_xillspec_table(tab, i0, i1, i2, i3, i4 + 1, i5);
  float* dat_1001 = get_xillspec_table(tab, i0, i1 + 1, i2, i3, i4 + 1, i5);
  float* dat_0101 = get_xillspec_table(tab, i0, i1, i2 + 1, i3, i4 + 1, i5);
  float* dat_0011 = get_xillspec_table(tab, i0, i1, i2, i3 + 1, i4 + 1, i5);

  float* dat_1101 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3, i4 + 1, i5);
  float* dat_1011 = get_xillspec_table(tab, i0, i1 + 1, i2, i3 + 1, i4 + 1, i5);
  float* dat_0111 = get_xillspec_table(tab, i0, i1, i2 + 1, i3 + 1, i4 + 1, i5);
  float* dat_1111 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4 + 1, i5);

  double f_0000 = (1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (1 - f4);
  double f_1000 = (f1) * (1.0 - f2) * (1.0 - f3) * (1 - f4);
  double f_0100 = (1.0 - f1) * (f2) * (1.0 - f3) * (1 - f4);
  double f_0010 = (1.0 - f1) * (1.0 - f2) * (f3) * (1 - f4);
  double f_0001 = (1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (f4);

  double f_1010 = (f1) * (1.0 - f2) * (f3) * (1 - f4);
  double f_1001 = (f1) * (1.0 - f2) * (1.0 - f3) * (f4);
  double f_0101 = (1.0 - f1) * (f2) * (1.0 - f3) * (f4);
  double f_1100 = (f1) * (f2) * (1.0 - f3) * (1 - f4);
  double f_0110 = (1.0 - f1) * (f2) * (f3) * (1 - f4);
  double f_0011 = (1.0 - f1) * (1.0 - f2) * (f3) * (f4);

  double f_1101 = (f1) * (f2) * (1.0 - f3) * (f4);
  double f_1110 = (f1) * (f2) * (f3) * (1 - f4);
  double f_1011 = (f1) * (1.0 - f2) * (f3) * (f4);
  double f_0111 = (1.0 - f1) * (f2) * (f3) * (f4);
  double f_1111 = (f1) * (f2) * (f3) * (f4);

  int ii;
  for (ii = 0; ii < n_ener; ii++)
  {
    flu[ii] =
      f_0000 * (double)dat_0000[ii] +
      f_1000 * (double)dat_1000[ii] +
      f_0100 * (double)dat_0100[ii] +
      f_0010 * (double)dat_0010[ii] +
      f_1100 * (double)dat_1100[ii] +
      f_1010 * (double)dat_1010[ii] +
      f_0110 * (double)dat_0110[ii] +
      f_1110 * (double)dat_1110[ii] +
      f_0001 * (double)dat_0001[ii] +
      f_1001 * (double)dat_1001[ii] +
      f_0101 * (double)dat_0101[ii] +
      f_0011 * (double)dat_0011[ii] +
      f_1101 * (double)dat_1101[ii] +
      f_1011 * (double)dat_1011[ii] +
      f_0111 * (double)dat_0111[ii] +
      f_1111 * (double)dat_1111[ii];
  }
}

static void interp_5d_tab(const xillTable* tab, double* flu, int n_ener,
                          double f1, double f2, double f3, double f4, double f5,
                          int i0, int i1, int i2, int i3, int i4,
                          int i5)
{
  float* dat_00000 = get_xillspec_table(tab, i0, i1, i2, i3, i4, i5);
  float* dat_10000 = get_xillspec_table(tab, i0, i1 + 1, i2, i3, i4, i5);
  float* dat_01000 = get_xillspec_table(tab, i0, i1, i2 + 1, i3, i4, i5);
  float* dat_00100 = get_xillspec_table(tab, i0, i1, i2, i3 + 1, i4, i5);

  float* dat_11000 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3, i4, i5);
  float* dat_10100 = get_xillspec_table(tab, i0, i1 + 1, i2, i3 + 1, i4, i5);
  float* dat_01100 = get_xillspec_table(tab, i0, i1, i2 + 1, i3 + 1, i4, i5);
  float* dat_11100 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4, i5);

  float* dat_00010 = get_xillspec_table(tab, i0, i1, i2, i3, i4 + 1, i5);
  float* dat_10010 = get_xillspec_table(tab, i0, i1 + 1, i2, i3, i4 + 1, i5);
  float* dat_01010 = get_xillspec_table(tab, i0, i1, i2 + 1, i3, i4 + 1, i5);
  float* dat_00110 = get_xillspec_table(tab, i0, i1, i2, i3 + 1, i4 + 1, i5);

  float* dat_11010 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3, i4 + 1, i5);
  float* dat_10110 = get_xillspec_table(tab, i0, i1 + 1, i2, i3 + 1, i4 + 1, i5);
  float* dat_01110 = get_xillspec_table(tab, i0, i1, i2 + 1, i3 + 1, i4 + 1, i5);
  float* dat_11110 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4 + 1, i5);

  float* dat_00001 = get_xillspec_table(tab, i0, i1, i2, i3, i4, i5 + 1);
  float* dat_10001 = get_xillspec_table(tab, i0, i1 + 1, i2, i3, i4, i5 + 1);
  float* dat_01001 = get_xillspec_table(tab, i0, i1, i2 + 1, i3, i4, i5 + 1);
  float* dat_00101 = get_xillspec_table(tab, i0, i1, i2, i3 + 1, i4, i5 + 1);

  float* dat_11001 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3, i4, i5 + 1);
  float* dat_10101 = get_xillspec_table(tab, i0, i1 + 1, i2, i3 + 1, i4, i5 + 1);
  float* dat_01101 = get_xillspec_table(tab, i0, i1, i2 + 1, i3 + 1, i4, i5 + 1);
  float* dat_11101 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4, i5 + 1);

  float* dat_00011 = get_xillspec_table(tab, i0, i1, i2, i3, i4 + 1, i5 + 1);
  float* dat_10011 = get_xillspec_table(tab, i0, i1 + 1, i2, i3, i4 + 1, i5 + 1);
  float* dat_01011 = get_xillspec_table(tab, i0, i1, i2 + 1, i3, i4 + 1, i5 + 1);
  float* dat_00111 = get_xillspec_table(tab, i0, i1, i2, i3 + 1, i4 + 1, i5 + 1);

  float* dat_11011 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3, i4 + 1, i5 + 1);
  float* dat_10111 = get_xillspec_table(tab, i0, i1 + 1, i2, i3 + 1, i4 + 1, i5 + 1);
  float* dat_01111 = get_xillspec_table(tab, i0, i1, i2 + 1, i3 + 1, i4 + 1, i5 + 1);
  float* dat_11111 = get_xillspec_table(tab, i0, i1 + 1, i2 + 1, i3 + 1, i4 + 1, i5 + 1);

  double f_00000 = (1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (1 - f4) * (1 - f5);
  double f_10000 = (f1) * (1.0 - f2) * (1.0 - f3) * (1 - f4) * (1 - f5);
  double f_01000 = (1.0 - f1) * (f2) * (1.0 - f3) * (1 - f4) * (1 - f5);
  double f_00100 = (1.0 - f1) * (1.0 - f2) * (f3) * (1 - f4) * (1 - f5);
  double f_00010 = (1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (f4) * (1 - f5);
  double f_10100 = (f1) * (1.0 - f2) * (f3) * (1 - f4) * (1 - f5);
  double f_10010 = (f1) * (1.0 - f2) * (1.0 - f3) * (f4) * (1 - f5);
  double f_01010 = (1.0 - f1) * (f2) * (1.0 - f3) * (f4) * (1 - f5);
  double f_11000 = (f1) * (f2) * (1.0 - f3) * (1 - f4) * (1 - f5);
  double f_01100 = (1.0 - f1) * (f2) * (f3) * (1 - f4) * (1 - f5);
  double f_00110 = (1.0 - f1) * (1.0 - f2) * (f3) * (f4) * (1 - f5);
  double f_11010 = (f1) * (f2) * (1.0 - f3) * (f4) * (1 - f5);
  double f_11100 = (f1) * (f2) * (f3) * (1 - f4) * (1 - f5);
  double f_10110 = (f1) * (1.0 - f2) * (f3) * (f4) * (1 - f5);
  double f_01110 = (1.0 - f1) * (f2) * (f3) * (f4) * (1 - f5);
  double f_11110 = (f1) * (f2) * (f3) * (f4) * (1 - f5);

  double f_00001 = (1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (1 - f4) * (f5);
  double f_10001 = (f1) * (1.0 - f2) * (1.0 - f3) * (1 - f4) * (f5);
  double f_01001 = (1.0 - f1) * (f2) * (1.0 - f3) * (1 - f4) * (f5);
  double f_00101 = (1.0 - f1) * (1.0 - f2) * (f3) * (1 - f4) * (f5);
  double f_00011 = (1.0 - f1) * (1.0 - f2) * (1.0 - f3) * (f4) * (f5);
  double f_10101 = (f1) * (1.0 - f2) * (f3) * (1 - f4) * (f5);
  double f_10011 = (f1) * (1.0 - f2) * (1.0 - f3) * (f4) * (f5);
  double f_01011 = (1.0 - f1) * (f2) * (1.0 - f3) * (f4) * (f5);
  double f_11001 = (f1) * (f2) * (1.0 - f3) * (1 - f4) * (f5);
  double f_01101 = (1.0 - f1) * (f2) * (f3) * (1 - f4) * (f5);
  double f_00111 = (1.0 - f1) * (1.0 - f2) * (f3) * (f4) * (f5);
  double f_11011 = (f1) * (f2) * (1.0 - f3) * (f4) * (f5);
  double f_11101 = (f1) * (f2) * (f3) * (1 - f4) * (f5);
  double f_10111 = (f1) * (1.0 - f2) * (f3) * (f4) * (f5);
  double f_01111 = (1.0 - f1) * (f2) * (f3) * (f4) * (f5);
  double f_11111 = (f1) * (f2) * (f3) * (f4) * (f5);

  int ii;
  for (ii = 0; ii < n_ener; ii++)
  {
    flu[ii] =
      f_00000 * (double)dat_00000[ii] +
      f_10000 * (double)dat_10000[ii] +
      f_01000 * (double)dat_01000[ii] +
      f_00100 * (double)dat_00100[ii] +
      f_11000 * (double)dat_11000[ii] +
      f_10100 * (double)dat_10100[ii] +
      f_01100 * (double)dat_01100[ii] +
      f_11100 * (double)dat_11100[ii] +
      f_00010 * (double)dat_00010[ii] +
      f_10010 * (double)dat_10010[ii] +
      f_01010 * (double)dat_01010[ii] +
      f_00110 * (double)dat_00110[ii] +
      f_11010 * (double)dat_11010[ii] +
      f_10110 * (double)dat_10110[ii] +
      f_01110 * (double)dat_01110[ii] +
      f_11110 * (double)dat_11110[ii]
      +
      f_00001 * (double)dat_00001[ii] +
      f_10001 * (double)dat_10001[ii] +
      f_01001 * (double)dat_01001[ii] +
      f_00101 * (double)dat_00101[ii] +
      f_11001 * (double)dat_11001[ii] +
      f_10101 * (double)dat_10101[ii] +
      f_01101 * (double)dat_01101[ii] +
      f_11101 * (double)dat_11101[ii] +
      f_00011 * (double)dat_00011[ii] +
      f_10011 * (double)dat_10011[ii] +
      f_01011 * (double)dat_01011[ii] +
      f_00111 * (double)dat_00111[ii] +
      f_11011 * (double)dat_11011[ii] +
      f_10111 * (double)dat_10111[ii] +
      f_01111 * (double)dat_01111[ii] +
      f_11111 * (double)dat_11111[ii];
  }
}

static void interp_2d_tab_incl(const xillTable* tab, double* flu, int n_ener,
                               const double* fac, int nfac, const std::vector<int>& ind, int iincl)
{
  assert(nfac == 2);
  assert(fac != nullptr);

  int i4 = ind[0];
  int i5 = ind[1];

  float* dat_0 = get_xillspec_table(tab, 0, 0, 0, 0, i4, i5);
  float* dat_1 = get_xillspec_table(tab, 0, 0, 0, 0, i4 + 1, i5);

  int ii;
  for (ii = 0; ii < n_ener; ii++)
  {
    flu[ii] = interp_lin_1d(fac[0], dat_0[ii], dat_1[ii]);
  }
}


static void interp_6d_tab_incl(const xillTable* tab, double* flu, int n_ener,
                               const double* fac, int nfac, const std::vector<int>& ind, int iincl)
{
  assert(nfac == 6);

  double s1[n_ener];
  double s2[n_ener];

  interp_5d_tab_incl(tab, s1, n_ener,
                     fac[1], fac[2], fac[3], fac[4], ind[0],
                     ind[1], ind[2], ind[3], ind[4], iincl);

  interp_5d_tab_incl(tab, s2, n_ener,
                     fac[1], fac[2], fac[3], fac[4], ind[0] + 1,
                     ind[1], ind[2], ind[3], ind[4], iincl);

  int ii;
  for (ii = 0; ii < n_ener; ii++)
  {
    flu[ii] = interp_lin_1d(fac[0], s1[ii], s2[ii]);
  }
}

static void interp_6d_tab(const xillTable* tab, double* flu, int n_ener,
                          const double* fac, int nfac, const std::vector<int>& ind)
{
  assert(nfac == 6);

  double s1[n_ener];
  double s2[n_ener];

  interp_5d_tab(tab, s1, n_ener,
                fac[1], fac[2], fac[3], fac[4], fac[5],
                ind[0],
                ind[1], ind[2], ind[3], ind[4], ind[5]);

  interp_5d_tab(tab, s2, n_ener,
                fac[1], fac[2], fac[3], fac[4], fac[5],
                ind[0] + 1,
                ind[1], ind[2], ind[3], ind[4], ind[5]);

  int ii;
  for (ii = 0; ii < n_ener; ii++)
  {
    flu[ii] = interp_lin_1d(fac[0], s1[ii], s2[ii]);
  }
}

/**
 * @brief check the boundary of the Ecut parameter and set ipol_factor accordingly
 * @detail grav. redshift can lead to the code asking for an Ecut not tabulated,
 * although actually observed ecut is larger. We set it to the lowest/highest value in this case
 * @param: (input) tab
 * @param: (input) param
 * @param: (output) ipol_fac
 **/
static void ensure_ecut_within_boundarys(const xillTable* tab, const xillTableParam* param, double* ipol_fac)
{
  // first need to make sure ECUT is a parameter in the table
  int index_ect = get_xilltab_param_index(tab, PARAM_ECT);
  if (index_ect >= 0)
  {
    if (param->ect <= tab->param_vals[tab->param_index[index_ect]][0])
    {
      ipol_fac[index_ect] = 0.0;
    }
    // can happen due to grav. redshift, although actually observed ecut is larger
    if (param->ect >=
      tab->param_vals[tab->param_index[index_ect]][tab->num_param_vals[tab->param_index[index_ect]] - 1])
    {
      ipol_fac[index_ect] = 1.0;
    }
  }
}

static void resetInpvalsToBoundaries(char* pname, float* p_inpVal, float tabValLo, float tabValHi)
{
  float inpVal = *p_inpVal;
  if (inpVal < tabValLo)
  {
    if (is_debug_run())
    {
      printf(" *** warning: paramter %s=%e below lowest table value, restetting to %e\n",
             pname, inpVal, tabValLo);
    }
    inpVal = tabValLo;
  }
  else if (inpVal > tabValHi)
  {
    if (is_debug_run())
    {
      printf("\n *** warning: paramter %s=%e above largest table value, restetting to %e\n",
             pname, inpVal, tabValHi);
    }
    inpVal = tabValHi;
  }

  *p_inpVal = inpVal;
}


xillSpec* interp_xill_table(const xillTable* tab, const xillTableParam* param, const std::vector<int> ind, int* status)
{
  CHECK_STATUS_RET(*status, nullptr);

  xillSpec* spec = nullptr;
  if (is_xill_model(param->model_type))
  {
    spec = new_xill_spec(1, tab->n_ener, status);
  }
  else
  {
    spec = new_xill_spec(tab->n_incl, tab->n_ener, status);
  }
  if (*status != EXIT_SUCCESS)
  {
    free_xill_spec(spec);
    return nullptr;
  }

  assert(spec != nullptr);
  assert(spec->n_ener == tab->n_ener);

  // set the energy grid
  int ii;
  for (ii = 0; ii < spec->n_ener; ii++)
  {
    spec->ener[ii] = tab->elo[ii];
  }
  spec->ener[spec->n_ener] = tab->ehi[spec->n_ener - 1];

  // set the inclination grid
  for (ii = 0; ii < spec->n_incl; ii++)
  {
    spec->incl[ii] = tab->incl[ii];
  }

  int nfac = tab->num_param;
  float* inp_param_vals = nullptr;
  double* ipol_fac = nullptr;
  do
  {
    inp_param_vals = get_xilltab_global_paramvals(param, status);
    CHECK_MALLOC_BREAK_STATUS(inp_param_vals, status)
    ipol_fac = (double*)malloc(sizeof(double) * nfac);
    CHECK_MALLOC_BREAK_STATUS(ipol_fac, status)
  }
  while (0);

  if (*status != EXIT_SUCCESS || inp_param_vals == nullptr || ipol_fac == nullptr)
  {
    free_xill_spec(spec);
    return nullptr;
  }


  /* calculate the interpolation factor for all parameters
   * ([nfac-1] is inclination, which might not be used) */
  int pind;
  for (ii = 0; ii < nfac; ii++)
  {
    // need the index
    pind = tab->param_index[ii];

    if (!((!is_xill_model(param->model_type)) && (tab->param_index[ii] == PARAM_INC)))
    {
      resetInpvalsToBoundaries(tab->param_names[ii], &(inp_param_vals[pind]), tab->param_vals[ii][0],
                               tab->param_vals[ii][tab->num_param_vals[ii] - 1]);
    }

    ipol_fac[ii] = (inp_param_vals[pind] - tab->param_vals[ii][ind[ii]]) /
      (tab->param_vals[ii][ind[ii] + 1] - tab->param_vals[ii][ind[ii]]);
  }

  free(inp_param_vals);

  // check the boundary of the Ecut parameter and set ipol_fac[ECUT] accordingly
  // (can happen due to strong grav. redshift)
  ensure_ecut_within_boundarys(tab, param, ipol_fac);

  if (tab->num_param == 5)
  {
    if (is_xill_model(param->model_type))
    {
      interp_5d_tab(tab, spec->flu[0], spec->n_ener,
                    ipol_fac[0], ipol_fac[1], ipol_fac[2], ipol_fac[3], ipol_fac[4], 0,
                    ind[0], ind[1], ind[2], ind[3], ind[4]);
    }
    else
    {
      // do not interpolate over the inclination (so skip the last parameter)
      nfac--;
      assert(nfac == 4);
      // get the spectrum for EACH flux bin
      for (ii = 0; ii < spec->n_incl; ii++)
      {
        interp_5d_tab_incl(tab, spec->flu[ii], spec->n_ener,
                           ipol_fac[0], ipol_fac[1], ipol_fac[2], ipol_fac[3], 0,
                           ind[0], ind[1], ind[2], ind[3], ii);
      }
    }
  }
  else if (tab->num_param == 6)
  {
    if (is_xill_model(param->model_type))
    {
      for (ii = 0; ii < spec->n_incl; ii++)
      {
        interp_6d_tab(tab, spec->flu[ii], spec->n_ener,
                      ipol_fac, nfac, ind);
      }
    }
    else
    {
      for (ii = 0; ii < spec->n_incl; ii++)
      {
        interp_6d_tab_incl(tab, spec->flu[ii], spec->n_ener,
                           ipol_fac, nfac, ind, ii);
      }
    }
  }
  else if (tab->num_param == 2)
  {
    if (is_xill_model(param->model_type))
    {
      RELXILL_ERROR("xillver type model with 2 parameters does not exist", status);
    }
    else
    {
      interp_2d_tab_incl(tab, spec->flu[ii], spec->n_ener,
                         ipol_fac, nfac, ind, ii);
    }
  }
  else
  {
    RELXILL_ERROR("Given number of parameters of the xillver table not implemented", status);
  }

  if (*status != EXIT_SUCCESS)
  {
    free_xill_spec(spec);
    return nullptr;
  }
  free(ipol_fac);

  return spec;
}


xillSpec *get_xillver_spectra_table(xillTableParam *param, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  xillTable *tab = nullptr;
  const char *fname = get_init_xillver_table(&tab, param->model_type, param->prim_type, status);
  CHECK_STATUS_AND_THROW(*status);

  assert(fname != nullptr);

  // =1=  get the inidices
  auto indparam = get_xilltab_indices_for_paramvals(param, tab, status);
  CHECK_STATUS_AND_THROW(*status);

  // =2=  check if the necessary spectra for interpolation are loaded
  check_xilltab_cache(fname, param, tab, indparam, status);
  CHECK_STATUS_AND_THROW(*status);

  // =3= interpolate values
  xillSpec *spec = interp_xill_table(tab, param, indparam, status);
  CHECK_STATUS_AND_THROW(*status);

  return spec;
}


/** @brief similar to get_xillver_spectra_table, but uses the full xill_param input (from xpsecc)
 *  @details see get_xillver_spectra_table for more details, this function is just a wrapper
 * @param param
 * @param status
 * @return xillSpec
 */
xillSpec* get_xillver_spectra(const xillParam* param, int* status)
{
  CHECK_STATUS_RET(*status, nullptr);
  xillTableParam* param_table = nullptr;
  xillSpec* xill_spec = nullptr;

  try
  {
    param_table = get_xilltab_param(param, status);
    CHECK_STATUS_AND_THROW(*status);
    xill_spec = get_xillver_spectra_table(param_table, status);
    free(param_table);
  }
  catch (const std::exception& e)
  {
    free(param_table);
    free_xill_spec(xill_spec);
    RELXILL_ERROR("getting spectra from xillver table failed", status);
    return nullptr;
  }

  return xill_spec;
}

/**
 * @brief calculates the angle averaged xillver spectrum from a xillSpec structure
 * @param xill_spec
 * @param status
 * @return (double*) spectrum
 */
double *calc_angle_averaged_xill_spec(const xillSpec *xill_spec, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  auto *output_spec = new double[xill_spec->n_ener];
  CHECK_MALLOC_RET_STATUS(output_spec, status, nullptr)

  memset(output_spec, 0, xill_spec->n_ener * sizeof(output_spec[0]));
  for (int ii = 0; ii < xill_spec->n_incl; ii++) {

    double incl_factor = norm_factor_semi_infinite_slab(xill_spec->incl[ii] * 180 / M_PI) / xill_spec->n_incl;
    for (int jj = 0; jj < xill_spec->n_ener; jj++) {
      output_spec[jj] += incl_factor * xill_spec->flu[ii][jj];
    }
  }

  return output_spec;
}

/**
 * @brief calculate the ENERGY flux in a given band from a given bin-integrated photon flux (cts/bin), the standard unit to store
 * all spectra in the relxill code
 **/
static double get_energy_flux_band(const double *ener,
                                   const double *photon_flux,
                                   int n_ener,
                                   double emin,
                                   double emax) {
  double sum = 0.0;
  for (int ii = 0; ii < n_ener; ii++) {
    if (ener[ii] >= emin && ener[ii + 1] <= emax) {
      sum += photon_flux[ii] * 0.5 * (ener[ii] + ener[ii + 1]);
    }
  }
  return sum;
}

/**
 * @brief calculate the PHOTON flux in a given band from a given bin-integrated photon flux (cts/bin), the standard unit to store
 * all spectra in the relxill code
 **/
static double get_photon_flux_band(const double *ener,
                                   const double *photon_flux,
                                   int n_ener,
                                   double emin,
                                   double emax) {
  double sum = 0.0;
  for (int ii = 0; ii < n_ener; ii++) {
    if (ener[ii] >= emin && ener[ii + 1] <= emax) {
      sum += photon_flux[ii];
    }
  }
  return sum;
}



EnerGrid *get_coarse_xillver_energrid() {

  if (global_xill_egrid_coarse == nullptr) {
    global_xill_egrid_coarse = new EnerGrid ; // new_EnerGrid(status);
    global_xill_egrid_coarse->nbins = N_ENER_COARSE;
    global_xill_egrid_coarse->ener = new double[global_xill_egrid_coarse->nbins + 1];
    get_log_grid(global_xill_egrid_coarse->ener, global_xill_egrid_coarse->nbins + 1,
                 EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION);
  }

  return global_xill_egrid_coarse;
}

double calcXillverNormFromPrimarySpectrum(const double *flux, const double *ener, const int n, int *status) {
  /* get the normalization of the spectrum with respect to xillver
 *  - everything is normalized using dens=10^15 cm^3
 *  - normalization defined, e.g., in Appendix of Dauser+2016
 *  - needs to be calculated on the specific energy grid (defined globally)
 */
  double keV2erg = 1.602177e-09;

  // need this to make sure no floating point problems arise
  double floatCompFactor = 1e-6;
  if (ener[n] < (EMAX_XILLVER_NORMALIZATION - floatCompFactor)
      || ener[0] > (EMIN_XILLVER_NORMALIZATION + floatCompFactor)) {
    RELXILL_ERROR("can not calculate the primary spectrum normalization", status);
    printf("  the given energy grid from %e keV to %e keV does not cover the boundaries\n", ener[0], ener[n]);
    printf("  from [%e,%e] necessary for the calculation\n", EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION);
    return 0.0;
  }

  double sum_pl = 0.0;
  for (int ii = 0; ii < n; ii++) {
    if (ener[ii] >= EMIN_XILLVER_NORMALIZATION && ener[ii] <= EMAX_XILLVER_NORMALIZATION) {
      sum_pl += flux[ii] * 0.5 * (ener[ii] + ener[ii + 1]) * 1e20 * keV2erg; // * 1e20 * kev... can be out of for loop to save a bit of computation
    }
  }
  double norm_xillver_table = 1e15 / 4.0 / M_PI;
  return sum_pl / norm_xillver_table;
}

/**
 * @brief cutoff pl model in xpsec units
 * @param pl_flux_xill [return] photon flux in counts/bin
 * @param ener
 * @param n_ener
 * @param xill_param
 * @param ener_shift: spectrum shifted in energy by this value, normalization is kept constant (i.e, norm not shifted by g^Gamma)
 */
void spec_cutoffpl(double *pl_flux_xill,
                   const double *ener,
                   int n_ener,
                   double ect,
                   double gamma,
                   double ener_shift = 1.0) {
  /** note that in case of the nthcomp model Ecut is in the frame of the primary source
      but for the bkn_powerlaw it is given in the observer frame */

  double ecut = ect * ener_shift; // Important: Ecut is given in the frame of the observer

  for (int ii = 0; ii < n_ener; ii++) {
    double en = 0.5 * (ener[ii] + ener[ii + 1]);
    pl_flux_xill[ii] = exp(1.0 / ecut) *
        pow(en, -gamma) *
        exp(-en / ecut) *
        (ener[ii + 1] - ener[ii]);
  }
}

/**
* @brief nthcomp model in xpsec units
* @param pl_flux_xill [return] photon flux in counts/bin
* @param ener
* @param n_ener
* @param xill_param
* @param ener_shift  spectrum shifted in energy by this value, normalization is kept constant! (i.e, norm not shifted by g^Gamma)
*/
void spec_nthcomp(double *pl_flux_xill,
                  const double *ener,
                  int n_ener,
                  const xillTableParam *xill_param,
                  double ener_shift = 1.0) {
  double nthcomp_param[5];
  double kTe = xill_param->ect; // Important: kTe is given in the frame of the source
  double z = 1 / ener_shift - 1; // convert energy shift to redshift
  get_nthcomp_param(nthcomp_param, xill_param->gam, kTe, z);
  c_donthcomp(ener, n_ener, nthcomp_param, pl_flux_xill);
}

/**
* @brief blackbody model in xpsec units
* @param pl_flux_xill [return] photon flux in counts/bin
* @param ener
* @param n_ener
* @param xill_param
*/
void spec_blackbody(double *pl_flux_xill, const double *ener, int n_ener, const xillTableParam *xill_param) {
  double en;
  for (int ii = 0; ii < n_ener; ii++) {
    en = 0.5 * (ener[ii] + ener[ii + 1]);
    pl_flux_xill[ii] = en * en / (pow(xill_param->kTbb, 4) * (exp(en / xill_param->kTbb) - 1));
    pl_flux_xill[ii] *= (ener[ii + 1] - ener[ii]);
  }
}

/**
 * @brief calculate the primary spectrum, which can be shifted by the factor ener_shift_source_obs
 * @param pl_flux_xill: photon flux in xspec units (cts/bin)
 * @param ener_shift_source_obs: energy shift applied the spectrum
 */
void calc_primary_spectrum(double *pl_flux_xill,
                           const double *ener,
                           int n_ener,
                           const xillTableParam *xill_param,
                           int *status,
                           double ener_shift_source_obs) {

  CHECK_STATUS_VOID(*status);
  assert(ener_shift_source_obs >= 0.0);

  // ecut and kTe are given in the frame of the source
  if (xill_param->prim_type == PRIM_SPEC_ECUT) {
    spec_cutoffpl(pl_flux_xill, ener, n_ener, xill_param->ect, xill_param->gam, ener_shift_source_obs);

  } else if (xill_param->prim_type == PRIM_SPEC_NTHCOMP) {
    spec_nthcomp(pl_flux_xill, ener, n_ener, xill_param, ener_shift_source_obs);

  } else if (xill_param->prim_type == PRIM_SPEC_BB) {
    spec_blackbody(pl_flux_xill, ener, n_ener, xill_param);

  } else {
    RELXILL_ERROR("trying to calculate a primary continuum for an undefined primary spectral type (should not happen!)",
                  status);
  }
}

double *calc_normalized_primary_spectrum(const double *ener, int n_ener,
                                         const relParam *rel_param, const xillTableParam *xill_param, int *status) {

  // need to create a specific energy grid for the primary component normalization to
  // fulfill the XILLVER NORM condition (Dauser+2016)
  EnerGrid *egrid = get_coarse_xillver_energrid();
  auto prim_spec_source = new double[egrid->nbins];
  calc_primary_spectrum(prim_spec_source, egrid->ener, egrid->nbins, xill_param, status);

  // calculate the normalization at the primary source (i.e., not shifted by the energy shift)
  double const
      norm_factor_prim_spec = 1. /
                              calcXillverNormFromPrimarySpectrum(prim_spec_source, egrid->ener, egrid->nbins, status);
  delete[] prim_spec_source;

  // if we have the LP geometry defined, the spectrum will be different between the primary source and the observer
  // due to energy shift: need to calculate the normalization for the source, but the spectrum at the observer
  double energy_shift = 1.0;
  if (rel_param != nullptr && rel_param->emis_type == EMIS_TYPE_LP) {
    energy_shift = energy_shift_source_obs(rel_param);
  }
  // output flux
  auto prim_spec_observer = new double[n_ener];
  calc_primary_spectrum(prim_spec_observer, ener, n_ener, xill_param, status, energy_shift);

  // take the normalization factor into account
  for (int ii = 0; ii < n_ener; ii++) {
    prim_spec_observer[ii] *= norm_factor_prim_spec;
  }

  return prim_spec_observer;
}

/**
 * @brief calculate the flux ration between the xillver reflection spectrum and its
 * input spectrum in the 0.1-1000 keV energy range
 * @param flu
 * @param ener
 * @param n_ener
 * @param xill_table_param
 * @param status
 * @return
 */
double get_xillver_fluxcorr(double *flu, const double *ener, int n_ener,
                            const xillTableParam *xill_table_param, int *status) {

  // use a coarser energy grid to evaluate the direct spectrum, as we are only interested
  // in the total flux in the band
  auto egrid_coarse = get_coarse_xillver_energrid();
  auto direct_spec =
      calc_normalized_primary_spectrum(egrid_coarse->ener, egrid_coarse->nbins,
                                       nullptr, xill_table_param, status);

  double const emin = 0.1;
  double const emax = 1000;
  double const ratio_refl_direct =
      get_energy_flux_band(ener, flu, n_ener, emin, emax) /
          get_energy_flux_band(egrid_coarse->ener, direct_spec, egrid_coarse->nbins, emin, emax);

  delete[] direct_spec;
  return ratio_refl_direct;
}

/**
 * @brief calculate the normalization factor to convert the initial xillver primary spectrum
 * (given by xill_param_0) under the given energy shift, by taking the fixed xillver incident
 * spectrum normalization in the 0.1-1000keV range into account (Dauser+16, A1). Note that the
 * normalization of the primary spectrum therefore scales as 1/norm_flux_cgs.
 * @param energy_shift_source_disk [energy shift from the source to the disk]
 * @param xill_param_0 [xillver parameters for the inital spectrum]
 * @return norm_flux_cgs
 */
double calc_xillver_normalization_change(double energy_shift, const xillTableParam *xill_param_0) {

  int status = EXIT_SUCCESS;
  EnerGrid *egrid = get_coarse_xillver_energrid();
  auto prime_spec = new double[egrid->nbins];

  calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_0, &status);
  double source_spec_norm_factor = 1. /
                                   calcXillverNormFromPrimarySpectrum(prime_spec, egrid->ener, egrid->nbins, &status);
  // normalized source spec: prime_source_source*source_spec_norm_factor

  // get parameters for the source spectrum (ecut shifted in energy)
  auto xill_param_source = new xillTableParam;
  *xill_param_source = *xill_param_0; // shallow copy
  xill_param_source->ect = xill_param_0->ect * energy_shift;

  // calculate the primary spectrum at the source
  calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_source, &status);
  double disk_spec_norm_factor =
      1. / calcXillverNormFromPrimarySpectrum(prime_spec, egrid->ener, egrid->nbins, &status);
  // normalized prime spec: prime_spec_disk*disk_spec_norm_factor

  delete xill_param_source;
  delete[] prime_spec;

  return disk_spec_norm_factor / source_spec_norm_factor;
}

/**
 * @brief calculate the normalization factor to convert the initial xillver primary spectrum
 * (given by xill_param_0) under the given energy shift, by taking the fixed xillver incident
 * spectrum normalization in the 0.1-1000keV range into account (Dauser+16, A1).
 * @param energy_shift_source_disk [energy shift from the source to the disk for each zone]
 * @param n_zones [number of zones]
 * @param xill_param_0 [xillver parameters for the inital spectrum]
 * @return norm_flux_cgs
 */
double *calc_xillver_normalization_change_source_to_disk(const double *energy_shift,
                                                         int n_zones,
                                                         const xillTableParam *xill_param_0) {

  int status = EXIT_SUCCESS;
  EnerGrid *egrid = get_coarse_xillver_energrid();
  auto prime_spec = new double[egrid->nbins];

  calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_0, &status);
  double source_spec_norm_factor = 1. /
                                   calcXillverNormFromPrimarySpectrum(prime_spec, egrid->ener, egrid->nbins, &status);
  // normalized disk spec: prime_spec_source*source_spec_norm_factor

  auto norm_factors = new double[n_zones];

  // get parameters for the source spectrum (ecut shifted in energy)
  auto xill_param_disk = new xillTableParam;
  *xill_param_disk = *xill_param_0; // shallow copy
  for (int ii = 0; ii < n_zones; ii++) {
    xill_param_disk->ect = xill_param_0->ect * energy_shift[ii];

    // calculate the primary spectrum at the source
    calc_primary_spectrum(prime_spec, egrid->ener, egrid->nbins, xill_param_disk, &status);
    double disk_spec_norm_factor = 1. /
                                   calcXillverNormFromPrimarySpectrum(prime_spec, egrid->ener, egrid->nbins, &status);
    // normalized prime spec: prime_spec_source*source_spec_norm_factor

    norm_factors[ii] = disk_spec_norm_factor / source_spec_norm_factor;
  }

  delete xill_param_disk;
  delete[] prime_spec;

  return norm_factors;
}

/**
 * @brief calculate the gshift flux correction factor
 * @detail The gshift-flux-correction factor is defined as the ratio of the flux-change of the xillver spectrum
 * by shifting it to a factor 1.5 lower energies wrt to the flux change of a powerlaw-spectrum (with the same
 * gamma), shifted by the same amoount.
 *
 * C_g = ( F(g=g_r)/F(g=g_0)  )  /  ( (g_r)^Gamma / (g_0)^Gamma )
 *     = ( F(g=1.5)/F(g=1) ) /  ( (1.5)^Gamma
 * @param  flux (double*): xillver spectrum
 * @param  energy (double*): energy grid on which the spectrum is defined
 * @param  n_energy (int)
 * @param  gamma (photon index)
 * @return fac_gshift_fluxcorr: ratio of xillver spectrum shifted to lower energy by a factor 1.5 (g=2/3) to a
 * shifted powerlaw with gamma (gamma is the input to the xillver spec)
 */
double get_xillver_gshift_fluxcorr(double *flu, const double *ener, int n_ener, double gamma) {

  const double gshift_refvalue = 2. / 3.;  // shift it 1.5 to lower energies

  auto ener_z = new double[n_ener + 1];
  for (int ii = 0; ii <= n_ener; ii++) {  // ener array has n_ener+1 entries
    ener_z[ii] = ener[ii] / gshift_refvalue;
  }

  auto flu_z = new double[n_ener];

  _rebin_spectrum(ener_z, flu_z, n_ener, ener, flu, n_ener);
  for (int ii = 0; ii < n_ener; ii++) {
    flu_z[ii] *= gshift_refvalue;  // take time dilation into account (dE already taken into account as bin-integ)
  }

  delete[] ener_z;

  double emin = 0.15;
  double emax = 500.0;

  // ratio of no-shift wrt to a shift if 1.5 to lower energies, i.e.
  // this is equivalent to F(g=1.5)/F(g=1)
  double gshift_ratio =
      get_photon_flux_band(ener, flu, n_ener, emin, emax) /
          get_photon_flux_band(ener, flu_z, n_ener, emin, emax);

  delete[] flu_z;

  return gshift_ratio / pow(1.5, gamma);
}

/**
 * @brief calculate the (energy) flux correction factors for a xillver spectrum.
 * @detail Two correction factors are calculated, which are both the ratio of the energy flux of two spectra:
 * (1) fac_fluxcorr: ratio of reflected xillver spectrum to its input spectrum (below 1 for low ionization and above 1
 * for high ionization)
 * (2) fac_gshift_fluxcorr: ratio of xillver spectrum shifted to lower energy by a factor 1.5 (g=2/3) to a
 * shifted powerlaw with gamma (gamma is the input to the xillver spec)
 * @param (output) fac_fluxcorr (double*)
 * @param (output) fac_gshift_fluxcorr (double*)
 * @param xill_param
 */
void get_xillver_fluxcorrection_factors(const xillSpec *xill_spec,
                                        double *fac_fluxcorr,
                                        double *fac_gshift_fluxcorr,
                                        xillTableParam *xill_table_param,
                                        int *status) {
  CHECK_STATUS_VOID(*status);

  double *angle_averaged_xill_spec = calc_angle_averaged_xill_spec(xill_spec, status);
  CHECK_STATUS_VOID(*status);

  if (fac_fluxcorr != nullptr) {
    *fac_fluxcorr =
        get_xillver_fluxcorr(angle_averaged_xill_spec, xill_spec->ener, xill_spec->n_ener, xill_table_param, status);
  }

  if (fac_gshift_fluxcorr != nullptr) {
    *fac_gshift_fluxcorr =
        get_xillver_gshift_fluxcorr(angle_averaged_xill_spec, xill_spec->ener, xill_spec->n_ener, xill_table_param->gam);
  }

  delete[] angle_averaged_xill_spec;

  if (*status != EXIT_SUCCESS) {
    RELXILL_ERROR("failed to calculate the flux correction factor", status);
  }

}

double norm_factor_semi_infinite_slab(double incl_deg) {
  return 0.5 * cos(incl_deg * M_PI / 180);
}

/**
 * @description: adds the proper flux normalization for a semi-infinate slab
 *  under inclination angle incl
 **/
void norm_xillver_spec(xillSpec *spec, double incl) {

  assert(spec->n_incl == 1); // make sure there is only one spectrum given (for the chosen inclination)

  double norm_factor = norm_factor_semi_infinite_slab(incl);
  for (int ii = 0; ii < spec->n_ener; ii++) {
    spec->flu[0][ii] *= norm_factor;
  }

}




/**
 * @Function: calc_xillver_angdep
 * @Synopsis: Calculate the Angle Weighted Xillver Spectrum on the Standard Relxill Spectrum
 * @Input: xill_spec
 *         rel_dist[n_incl]
 * @Output: xill_flux  (needs to be allocated, will be overwritten)
 **/
void calc_xillver_angdep(double *xill_flux, xillSpec *xill_spec, const double *dist, const int *status) {

  CHECK_STATUS_VOID(*status);

  int ii;
  int jj;
  for (ii = 0; ii < xill_spec->n_ener; ii++) {
    xill_flux[ii] = 0.0;
  }

  for (ii = 0; ii < xill_spec->n_incl; ii++) {
    for (jj = 0; jj < xill_spec->n_ener; jj++) {
      xill_flux[jj] += dist[ii] * xill_spec->flu[ii][jj];
    }
  }

}


/**
 * @Function: calculateAngleWeightedXillverSpectrum
 * @Synopsis: Calculate the Angle Weighted Xillver Spectrum on the Standard Relxill Spectrum
 * @Input:  ener[n_ener]
 *         rel_dist[n_incl]
 *         xill_spec
 * @Output: o_xill_flux  (needs to be allocated, will be overwritten)
 *  [reason for the required allocation is that this will be called in a large
 *   loop and otherwise we would need to allocate a 3000 bin array very often]
 **/
void calculateAngleWeightedXillverSpectrum(double *o_xill_flux,
                                           double *ener,
                                           int n_ener,
                                           double *rel_dist,
                                           xillSpec *xill_spec,
                                           int *status) {
  CHECK_STATUS_VOID(*status);

  auto xill_angdist_inp = new double[xill_spec->n_ener];

  calc_xillver_angdep(xill_angdist_inp, xill_spec, rel_dist, status);

  _rebin_spectrum(ener, o_xill_flux, n_ener,
                  xill_spec->ener, xill_angdist_inp, xill_spec->n_ener);

  delete[] xill_angdist_inp;

}



/**
 * calculate a xillver spectrum and combine it according to the angle distribution
 * (stored in rel_cosne_dist)
 * @param [output] xill_flux  (needs to be allocated)
 * @param ener
 * @param n_ener
 * @param xill_param
 * @param rel_cosne_dist
 * @param status
 */
void getAnglecorrXillSpec(double *xill_flux, double *ener, int n_ener, xillParam *xill_param,
                          double *rel_cosne_dist, int *status) {

  CHECK_STATUS_VOID(*status);

  xillSpec *xillSpecAngles = get_xillver_spectra(xill_param, status);
  calculateAngleWeightedXillverSpectrum(xill_flux, ener, n_ener, rel_cosne_dist, xillSpecAngles, status);

}
