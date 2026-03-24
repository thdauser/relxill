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

    Copyright 2025 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef ENVOPTIONS_H_
#define ENVOPTIONS_H_

#include "common.h"

double get_env_otherwise_default(const char *env, double def_value);

/** get the relxill table path (dynamically from env variable)  **/
char *get_relxill_table_path(void);

/** check if we are currently debugging the model **/
int is_debug_run(void);

int shouldAuxInfoGetPrinted(void);

int shouldOutfilesBeWritten(void);

int get_from_env_rrad_bbody_order(void);

int get_from_env_rrad_bbody_tprofile(void);

int constantDiskDensity(void);

/** check if we should return the relline/relconv physical norm from ENV **/
int do_not_normalize_relline(void);

int get_from_env_relline_num_zones(void);

void print_env_options(void);

#endif /* ENVOPTIONS_H_ */