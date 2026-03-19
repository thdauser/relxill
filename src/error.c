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

#include "error.h"
#include <fitsio.h>

void relxill_error(const char *const func, const char *const msg, int *status) {
  *status = EXIT_FAILURE;
  printf(" *** error in relxill (%s): %s!\n", func, msg);
}

void relxill_warning(const char *const msg) {
  printf(" *** warning from relxill: %s!\n", msg);
}

/** check and report FITS error   */
void relxill_check_fits_error(const int *status) {
  if (*status != EXIT_SUCCESS) {
    char errtext[30];
    fits_get_errstatus(*status, errtext);
    printf("   cfitsio error: %s \n", errtext);
  }
}

void check_relxill_error(const char *const func, const char *const msg, int *status) {
  if (*status != EXIT_SUCCESS) {
    *status = EXIT_FAILURE;
    printf(" *** error in relxill (%s): %s!\n", func, msg);
  }
}

void print_relxill_test_msg(const char *const func, const char *const msg) {
  printf(" TEST: %s %s", func, msg);
}