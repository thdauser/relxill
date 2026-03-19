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

#include "Relxill_errorhandling.h"
#include "Xilltable_handling.h"
#include "Rebin.h"

#include "Types.h"

extern "C" {
#include "xilltable_reading.h"
#include "common.h"
#include "relutility.h"
}

#include <vector>
#include <memory>
#include <cassert>

/** @brief main routine for the xillver table: returns a spectrum for the given parameters
 *  @details
 *   - decides if the table needs to be initialized and/or more data loaded
 *   - automatically normalizes  the spectra to logN=1e15 and logXi=0
 *   - for relxill type models it returns spectra for all available inclinations
 * @param n_incl, n_ener
 * @param status
 * @return xillSpec
 */


xillSpec* new_xill_spec(int n_incl, int n_ener, int* status)
{
  CHECK_STATUS_RET(*status, nullptr);

  auto* spec = new xillSpec;
  CHECK_MALLOC_RET_STATUS(spec, status, nullptr)

  spec->n_ener = n_ener;
  spec->n_incl = n_incl;
  spec->ener = nullptr;
  spec->incl = nullptr;
  spec->flu = nullptr;

  try
  {
    spec->ener = new double[n_ener + 1];
    spec->incl = new double[n_incl];
    spec->flu = new double*[n_incl];

    for (int ii = 0; ii < n_incl; ii++)
    {
      spec->flu[ii] = new double[n_ener];
    }
  }  catch (const std::bad_alloc&) {
    *status = EXIT_FAILURE;
  }

  if (*status == EXIT_FAILURE)
  {
    printf(" error: allocating memeory for the xillver spectrum failed!\n");
    free_xill_spec(spec);
    return nullptr;
  }

  return spec;
}

void free_xill_spec(xillSpec* spec)
{
  if (spec != nullptr)
  {
    delete[] spec->ener;
    delete[] spec->incl;
    if (spec->flu != nullptr)
    {
      for (int ii = 0; ii < spec->n_incl; ii++)
      {
        delete[] spec->flu[ii];
      }
      delete[] spec->flu;
    }
    delete spec;
  }
}

/**
 * @brief Convert xillver table parameters to a global parameter values array
 * @details Creates an array indexed by global parameter identifiers (PARAM_GAM, PARAM_AFE, etc.)
 *          containing the corresponding values from the xillTableParam structure. This allows
 *          unified access to all xillver parameters using their global index constants.
 * @param param Input structure containing xillver table parameters
 * @return Pointer to the allocated float array with N_PARAM_MAX elements, or NULL on error.
 *         Caller is responsible for freeing the returned array.
 */
std::vector<float> get_xilltab_global_paramvals(const xillTableParam* param)
{
  std::vector<float> global_param_vals(N_PARAM_MAX);

  global_param_vals[PARAM_GAM] = static_cast<float>(param->gam);
  global_param_vals[PARAM_AFE] = static_cast<float>(param->afe);
  global_param_vals[PARAM_LXI] = static_cast<float>(param->lxi);
  global_param_vals[PARAM_ECT] = static_cast<float>(param->ect);
  global_param_vals[PARAM_DNS] = static_cast<float>(param->dens);
  global_param_vals[PARAM_KTB] = static_cast<float>(param->kTbb);
  global_param_vals[PARAM_FRA] = static_cast<float>(param->frac_pl_bb);
  global_param_vals[PARAM_RID] = static_cast<float>(param->radindex);
  global_param_vals[PARAM_INC] = static_cast<float>(param->incl);

  return global_param_vals;
}

std::vector<int> get_xilltab_indices_for_paramvals(const xillTableParam* param, const xillTable* tab, int* status)
{
  CHECK_STATUS_RET(*status, std::vector<int>{});

  // Use unique_ptr for automatic memory management
  auto inp_global_paramvals =
    get_xilltab_global_paramvals(param);

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
  VecI loop_start(max_dims);
  VecI loop_end(max_dims);

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
                               const std::vector<double>& fac, int ind_radzone, int ind_incl)
{
  // make sure this is only for 2D interpolation
  int dim = 2;
  assert(fac.size() == dim);
  assert(tab->num_param == dim);

  int i4 = ind_radzone;
  int i5 = ind_incl;

  float* dat_0 = get_xillspec_table(tab, 0, 0, 0, 0, i4, i5);
  float* dat_1 = get_xillspec_table(tab, 0, 0, 0, 0, i4 + 1, i5);

  int ii;
  for (ii = 0; ii < n_ener; ii++)
  {
    flu[ii] = interp_lin_1d(fac[0], dat_0[ii], dat_1[ii]);
  }
}


static void interp_6d_tab_incl(const xillTable* tab, double* flu, int n_ener,
                               const std::vector<double>& fac, const std::vector<int>& ind, int iincl)
{
  // make sure this is only for 2D interpolation
  int dim = 6;
  assert(fac.size() == dim);
  assert(ind.size() == dim);
  assert(tab->num_param == dim);

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
                          const std::vector<double>& fac, const std::vector<int>& ind)
{
  // make sure this is only for 2D interpolation
  int dim = 6;
  assert(fac.size() == dim);

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
static void ensure_ecut_within_boundarys(const xillTable* tab, const xillTableParam* param, std::vector<double> ipol_fac)
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


xillSpec* interp_xill_table(const xillTable* tab, const xillTableParam* param, const std::vector<int>& ind, int* status)
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
  auto inp_param_vals = get_xilltab_global_paramvals(param);
  std::vector<double> ipol_fac(nfac);

  if (*status != EXIT_SUCCESS)
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
                      ipol_fac, ind);
      }
    }
    else
    {
      for (ii = 0; ii < spec->n_incl; ii++)
      {
        interp_6d_tab_incl(tab, spec->flu[ii], spec->n_ener,
                           ipol_fac, ind, ii);
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
      for (ii = 0; ii < spec->n_incl; ii++)
      {
        interp_2d_tab_incl(tab, spec->flu[ii], spec->n_ener,
                           ipol_fac, ind[0], ii);
        spec->flu[ii][0] = 1.0;
      }
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

  // =2= check if the necessary spectra for interpolation are loaded
  check_xilltab_cache(fname, param, tab, indparam, status);
  CHECK_STATUS_AND_THROW(*status);

  // =3= interpolate values
  xillSpec *spec = interp_xill_table(tab, param, indparam, status);
  CHECK_STATUS_AND_THROW(*status);

  return spec;
}
