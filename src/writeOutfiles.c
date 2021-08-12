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
#include <stdio.h>
#include "writeOutfiles.h"

static void fclose_errormsg(FILE *fp, char *foutName) {
  if (fclose(fp)) {
    printf(" *** error : failed writing file %s \n", foutName);
  }

}

/**
 *  @function: write_binned_data_to_file
 *  @synopsis: write data to file, where the x-value is interpreted as a binned grid,
 *    with rlo[ii]=rad[ii] and rhi[ii]=rad[ii+1]
 *  @input:
 *   - double[n_rad+1] rad
 *   - double[n_rad] intens
 *   - int n_rad
 **/
void write_binned_data_to_file(char *foutName, const double *rad, double *intens, int n_rad) {
  FILE *fp = fopen(foutName, "w+");
  for (int ii = 0; ii < n_rad; ii++) {
    fprintf(fp, " %e \t %e \t %e \n", rad[ii], rad[ii + 1], intens[ii]);
  }
  fclose_errormsg(fp, foutName);
}

/**
 *  @function: write_data_to_file
 *  @synopsis: write data to file, where the x-value and y-value is written
 *  @input:
 *   - double[n_rad] rad
 *   - double[n_rad] intens
 *   - int n_rad
 **/
void write_data_to_file(const char *foutName, double *rad, double *intens, int n_rad) {

  FILE *fp = fopen(foutName, "w+");
  for (int ii = 0; ii < n_rad; ii++) {
    fprintf(fp, " %e \t %e \n", rad[ii], intens[ii]);
  }
  if (fclose(fp)) {
    exit(1);
  }
}

void save_relline_radial_flux_profile(double *rad, double *intens, int n_rad) {
  char *fname = "test_relline_radialFluxProfile.dat";
  assert(intens != NULL);
  write_data_to_file(fname, rad, intens, n_rad);
}

void save_relline_profile(relline_spec_multizone *spec) {
  char *fname = "test_relline_profile.dat";
  assert(spec != NULL);
  write_binned_data_to_file(fname, spec->ener, spec->flux[0], spec->n_ener);
}

void save_xillver_spectrum(const double *ener, double *flu, int n_ener, char *fname) {
  write_binned_data_to_file(fname, ener, flu, n_ener);
}

/**
 * @function save_emis_profiles
 * @param sysPar
 * @synopsis: saves the emissivity profile, the incident one and (if exists) also
 * the returning radiation profile
 */
void save_emis_profiles(RelSysPar *sysPar) {
  assert(sysPar->emis != NULL);
  write_data_to_file("test_emis_profile.dat",
                     sysPar->emis->re, sysPar->emis->emis, sysPar->emis->nr);
}

