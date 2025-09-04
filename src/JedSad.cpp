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

    Copyright 2024 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "JedSad.h"
// #include "mlinterp.hpp" // todo: check if this would be a faster way of interpolation

extern "C" {
#include "xilltable.h"
#include "relutility.h"
}

#define CHECK_STATUS_CFITSIO(status) \
 if (EXIT_SUCCESS!=status) {printf(" error in cfitsio function with code %i\n",status); throw std::exception();}



void JedsadTable::read_jstable_params(fitsfile* fptr){

  int extver = 0;
  string extname = "PARAMETERS";
  int status = EXIT_SUCCESS;
  fits_movenam_hdu_cpp(fptr, BINARY_TBL, extname.c_str(), extver, &status);
  CHECK_STATUS_CFITSIO(status)

  // we know the column numbers
  int colnum_n = 2;
  int colnum_vals = 3;

  long n;
  fits_get_num_rows(fptr, &n, &status);
  CHECK_STATUS_CFITSIO(status)

  if (num_params() != n) {
    printf("wrong format of the JEDSAD table (not the correct number of parameters)");
    throw std::exception();
  }

  int anynul = 0;
  double nullval = 0.0;

  int ii;
  char strnull[10];
  strcpy(strnull, " ");

  fits_read_col(fptr, TINT, colnum_n, 1, 1, n, &nullval, num_param_vals, &anynul, &status);
  CHECK_STATUS_CFITSIO(status)

  assert(param_vals.size() == n);
  for (ii = 0; ii < n; ii++) {
    /** here we load the parameter values **/
    param_vals[ii] = new double[num_param_vals[ii]];
    fits_read_col(fptr, TDOUBLE, colnum_vals, ii + 1, 1, num_param_vals[ii], &nullval, param_vals[ii], &anynul,&status);
    CHECK_STATUS_CFITSIO(status)
  }

  int n_data_entries = 1;
  for (ii = 0; ii < n; ii++){
    n_data_entries *= num_param_vals[ii];
  }

  // read the data
  string extname_data = "DATA";
  fits_movenam_hdu_cpp(fptr, BINARY_TBL, extname_data.c_str(), extver, &status);
  CHECK_STATUS_CFITSIO(status)

  int n_columns_data;
  fits_get_num_cols(fptr, &n_columns_data, &status);
  CHECK_STATUS_CFITSIO(status)
  assert(n_columns_data == raw_data.size());

  long n_rows_data;
  fits_get_num_rows(fptr, &n_rows_data, &status);
  CHECK_STATUS_CFITSIO(status)
  assert(n_rows_data == n_data_entries);
  m_num_data_vals = n_data_entries;

  for(ii = 0; ii< n_columns_data; ii++) {
    raw_data[ii] = new double[n_rows_data];
    fits_read_col(fptr, TDOUBLE, ii + 1, 1, 1, n_rows_data, &nullval, raw_data[ii], &anynul, &status);
  }

}


// from the index for each of the parameters, calculate the index of the data
int JedsadTable::get_data_index(std::vector<int> pind) {


  assert (num_params() == 9);

  return (((((((pind[0]
                * num_param_vals[1] + pind[1])
               * num_param_vals[2] + pind[2]) * num_param_vals[3] + pind[3])
             * num_param_vals[4] + pind[4]) * num_param_vals[5] + pind[5])
           * num_param_vals[6] + pind[6]) * num_param_vals[7] + pind[7])
         * num_param_vals[8] + pind[8];
}




void JedsadTable::read_table() {

  if (!raw_data.empty()) {
    printf(" data not empty, skipping reading the table\n");
    return;
  }

  // init sizes
  param_vals.resize(num_params());
  raw_data.resize(JedsadTableInformation::instance().n_data);

  int status = EXIT_SUCCESS;
  fitsfile *fptr = open_fits_table_stdpath(m_fullfilename.c_str(), &status);
  CHECK_STATUS_CFITSIO(status)

  assert(fptr != nullptr);
  read_jstable_params(fptr);

  fits_close_file(fptr, &status);
  CHECK_STATUS_CFITSIO(status)

}

std::vector<double> JedsadTable::interpolate(const double *param_array) {

  if (raw_data.empty()) {
    read_table();
  }

  const int n_params = num_params();
  const int *_num_param_values = num_param_vals;
  const std::vector<double *> _p_vals = param_vals;

  // Arrays to store interpolation indices and weights
  std::vector<int> lower_indices(n_params);
  std::vector<int> upper_indices(n_params);
  std::vector<double> weights(n_params);

  // Find bracketing indices and calculate weights for each parameter
  for (int dim = 0; dim < n_params; dim++) {
    double query_val = param_array[dim];
    double *param_grid = _p_vals[dim];
    int n_grid_points = _num_param_values[dim];

    // Use binary search to find the bracketing indices
    int found_index = binary_search(param_grid, n_grid_points, query_val);

    // Handle boundary cases and set interpolation indices
    if (query_val <= param_grid[0]) {
      // Below grid: use first point
      lower_indices[dim] = 0;
      upper_indices[dim] = 0;
      weights[dim] = 1.0;
    } else if (query_val >= param_grid[n_grid_points - 1]) {
      // Above grid: use last point
      lower_indices[dim] = n_grid_points - 1;
      upper_indices[dim] = n_grid_points - 1;
      weights[dim] = 1.0;
    } else if (found_index >= 0) {
      // Check if it's truly an exact match
      if (param_grid[found_index] == query_val) {
        // Exact match found
        lower_indices[dim] = found_index;
        upper_indices[dim] = found_index;
        weights[dim] = 1.0;
      } else {
        // Not an exact match, interpolate between found_index and found_index+1
        lower_indices[dim] = found_index;
        upper_indices[dim] = found_index + 1;

        // Calculate interpolation weight
        double lower_val = param_grid[lower_indices[dim]];
        double upper_val = param_grid[upper_indices[dim]];
        weights[dim] = (query_val - lower_val) / (upper_val - lower_val);

      }
    } else {
      // Not found: binary_search returns -(insertion_point + 1)
      // insertion_point is where the value would be inserted
      int insertion_point = -(found_index + 1);
      lower_indices[dim] = insertion_point - 1;
      upper_indices[dim] = insertion_point;

      // Ensure indices are within bounds
      if (lower_indices[dim] < 0) {
        lower_indices[dim] = 0;
        upper_indices[dim] = 0;
        weights[dim] = 1.0;
      } else if (upper_indices[dim] >= n_grid_points) {
        lower_indices[dim] = n_grid_points - 1;
        upper_indices[dim] = n_grid_points - 1;
        weights[dim] = 1.0;
      } else {
        // Calculate interpolation weight
        double lower_val = param_grid[lower_indices[dim]];
        double upper_val = param_grid[upper_indices[dim]];
        weights[dim] = (query_val - lower_val) / (upper_val - lower_val);
      }
    }
  }

  // Perform 9D multilinear interpolation
  // We need to iterate through all 2^9 = 512 corner points of the hypercube
  const int n_corners = 1 << n_params; // 2^9 = 512
  const int n_data_columns = raw_data.size();

  // Storage for interpolated results (one for each data column)
  std::vector<double> interpolated_values(n_data_columns, 0.0);

  for (int corner = 0; corner < n_corners; corner++) {
    double corner_weight = 1.0;
    std::vector<int> corner_indices(n_params);

    // Determine which corner of the hypercube we're at
    for (int dim = 0; dim < n_params; dim++) {
      bool use_upper = (corner >> dim) & 1;
      corner_indices[dim] = use_upper ? upper_indices[dim] : lower_indices[dim];
      double dim_weight = use_upper ? weights[dim] : (1.0 - weights[dim]);
      corner_weight *= dim_weight;
    }

    // Calculate flat index for this corner point
    int flat_index = get_data_index(corner_indices);

    // Add contribution from this corner to each data column
    for (int col = 0; col < n_data_columns; col++) {
      interpolated_values[col] += corner_weight * raw_data[col][flat_index];
    }
  }

  return interpolated_values;
}

PrimespecParams convert_jedsad_to_primaryspec_params(const ModelDefinition& model_definition){

  auto jedsadParams = JedsadParams(model_definition);

  double* param_array = jedsadParams.get_param_array();
  auto result = JedsadTable::instance().interpolate(param_array);
  delete[](param_array);


  return PrimespecParams(result[0], result[1]);
}


 // jedsadTable *cached_jedsad_table = nullptr;