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

    Copyright 2016 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELUTILITY_H_
#define RELUTILITY_H_

#include "relbase.h"


/****** DEFINE FUNCTION DEFINITIONS ******/

#define RELXILL_ERROR(msg,status) (relxill_error(__func__, msg,status))

#define CHECK_RELXILL_ERROR(msg,status) (check_relxill_error(__func__, msg,status))

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);

#define CHECK_STATUS_VOID(status) \
  if (EXIT_SUCCESS!=status) return;

#define CHECK_STATUS_BREAK(status) \
  if (EXIT_SUCCESS!=status) break;

#define CHECK_MALLOC_VOID_STATUS(a,status) \
	if (NULL==a) { \
		RELXILL_ERROR("memory allocation failed",status); \
		return;\
	}

#define CHECK_MALLOC_RET_STATUS(a,status,retval) \
	if (NULL==a) { \
		RELXILL_ERROR("memory allocation failed",status); \
		return retval;\
	}

#define CHECK_STATUS_RET(status, retval) \
  if (EXIT_SUCCESS!=status) return(retval);



/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double* ener, int n_ener, double emin, double emax);

/* get the current version number */
void get_version_number(char** vstr, int* status);

/* print relxill error message */
void relxill_error(const char* const func, const char* const msg, int* status);

/* check and print relxill error message */
void check_relxill_error(const char* const func, const char* const msg, int* status);

#endif /* RELUTILITY_H_ */
