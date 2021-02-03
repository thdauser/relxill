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

    Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "LocalModel.h"
#include "XspecSpectrum.h"
#include "xspec_wrapper_lmodels.h"

#include <vector>

extern "C" {
#include "relutility.h"
}



void default_eval_local_model(ModelName model_name){

  DefaultSpec default_spec{};
  LocalModel local_model(model_name);

  XspecSpectrum spec = default_spec.get_xspec_spectrum();
  local_model.eval_model(spec);

  assert(calcSum(default_spec.flux, default_spec.num_flux_bins) > 1e-6);

}


ModelName find_name_string_in_xspec_database(const XspecModelDatabase& database, const std::string& model_string){

  auto models = database.all_models();

  for (auto & model : models){
    if (model.second.name()==model_string){
      return model.first;
    }

  }
  throw ModelEvalFailed("failed to find a model with this name");
}


ModelName get_model_name_from_string(const std::string& model_string){

  XspecModelDatabase database{};
  std::cout << "  evaluating model: " << model_string << std::endl;
  return find_name_string_in_xspec_database(database, model_string);

}

int main(int argc, char *argv[]) {

    if (argc != 2)  {
      printf(" Missing argument: ");
      printf("  - Return version number: ./test_sta version");
      printf("  - Evaluate model: ./test_sta <model_name>");

    } else if (strcmp(argv[1],"version") == 0) {
      printf("%s",PROJECT_VER);

    } else {

      ModelName model_name = get_model_name_from_string(std::string(argv[1]) );
      default_eval_local_model(model_name);

    }

    return EXIT_SUCCESS;
}
