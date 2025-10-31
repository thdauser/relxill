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


#ifndef RELXILL_TYPES_H
#define RELXILL_TYPES_H

#include <memory>
#include <vector>
using VecD = std::vector<double>;
using VecD_2D = std::vector<std::vector<double>>;
using VecU_D = std::vector<std::unique_ptr<double[]>>;

using VecI = std::vector<int>;

#endif //RELXILL_TYPES_H