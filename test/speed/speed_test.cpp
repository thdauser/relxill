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

#include <chrono>

extern "C" {
#include "relutility.h"
}



void eval_local_model_param_range(ModelName model_name, XPar param, double pmin, double pmax, int npar){

  DefaultSpec default_spec{};
  XspecSpectrum spec = default_spec.get_xspec_spectrum();


  LocalModel local_model(model_name);

  for (int ii=0; ii<npar; ii++){
    double param_value = (pmax - pmin) * ( static_cast<double>(ii) / static_cast<double>(npar) );
    local_model.set_par(param, param_value);
    local_model.eval_model(spec);
  }

}


void eval_model_relat_param_changes(ModelName model_name, const int num_evaluations){
  XPar rel_param = XPar::a;
  eval_local_model_param_range(model_name, rel_param, 0.5, 0.998, num_evaluations);
}

void eval_model_xillver_param_changes(ModelName model_name, const int num_evaluations){
  XPar rel_param = XPar::logxi;
  // only choose a small difference here such that the table does not need to be re-loaded
  eval_local_model_param_range(model_name, rel_param, 3.1, 3.2, num_evaluations);
}


// ------------------------- //
int main(int argc, char *argv[]) {

  if (argc < 2)  {
    printf(" Missing argument: ");
    printf("./speed_test <model> [<rel|m_cache_xill>]");
  } else {

    const int num_evaluations = 100;

    const int num_zones = 25;
    setenv("RELXILL_NUM_RZONES", std::to_string(num_zones).c_str(), 1);
    printf(" testing %s for number of zones %i \n",argv[1], num_zones);

    ModelName model_name = ModelDatabase::instance().model_name(std::string(argv[1]));

    auto tstart = std::chrono::steady_clock::now();

    if (argc==2){
      assert(num_evaluations>1);
      eval_model_relat_param_changes(model_name, num_evaluations/2);
      eval_model_xillver_param_changes(model_name, num_evaluations/2);
    } else {

      if  (std::string(argv[2]) == "rel"){
        eval_model_relat_param_changes(model_name, num_evaluations);
      } else if (std::string(argv[2]) == "m_cache_xill") {
        eval_model_xillver_param_changes(model_name, num_evaluations);
      } else {
        std::cout<< " speed test error: given argument " << std::string(argv[2]) << " not known" << std::endl;
      }
    }


    auto time_elapsed_msec = std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::steady_clock::now() - tstart ).count();

    printf("Time elapsed:              %.2fsec\n", static_cast<double>(time_elapsed_msec) * 0.001);
    printf("Time per model evaluation: %.0fmsec ", static_cast<double>(time_elapsed_msec) / num_evaluations);

  }

  return EXIT_SUCCESS;
}
