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

#ifndef RELXILL__ERROR_H_
#define RELXILL__ERROR_H_

#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* print relxill error message */
void relxill_error(const char* func, const char* msg, int* status);

/* print a standardized warning message */
void relxill_warning(const char* msg);

/* check and print relxill error message */
void check_relxill_error(const char* func, const char* msg, int* status);

/** check and report FITS error   */
void relxill_check_fits_error(const int* status);

void print_relxill_test_msg(const char* func, const char* msg);

#ifdef __cplusplus
}
#endif

#endif //RELXILL__ERROR_H_