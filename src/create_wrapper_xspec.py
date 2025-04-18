#! /usr/bin/env python3
import re
import sys

# Tool to Automatically Construct the Cpp-Model-Wrapper from the Xspec lmodel.dat file
# - models are separated by an empty line
# - for <model_name> the C-function has to be called c_<local_model_prefix><model_name>

local_model_prefix = "lmod"

header = """/*
   *** AUTOMATICALLY GENERATED: DO NOT EDIT THIS FILE ***   

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
       
    *** AUTOMATICALLY GENERATED: DO NOT EDIT THIS FILE ***       
*/

"""


class ModelDefinition:
    def __init__(self, model_definition_string):
        self.function_name = None
        self.model_name = None
        self.params = None

        self.set_model_and_function_name(model_definition_string)
        self.params = self.get_param_list(model_definition_string)

        self.test_definition()

    def get_param_list(self, definition_string):
        if self.model_name is not None:
            return parse_param_list(definition_string)
        else:
            print(" *** error ***: could not parse the following model definition")
            print(definition_string)
            return None

    def set_model_and_function_name(self, definition_string):
        first_line = split_in_lines_and_remove_empty(definition_string)[0]

        res = re.match(rf'.*c_({local_model_prefix}\w+).*', first_line)
        if res is not None:
            self.function_name = res.groups()[0]
            self.model_name = first_line.split(' ')[0]

    def test_definition(self):
        if self.params is not None:
            print("    - " + self.model_name + "  found " + str(len(self.params)) +
                  " parameters \t  (function name: " + self.function_name + ")")


def remove_empty_strings(list_strings):
    for string in list_strings:
        if len(string) <= 1:
            list_strings.remove(string)

    return list_strings


def remove_empty_line(list_strings):
    for string in list_strings:
        string.lstrip()
        print(string)

    return list_strings


def read_file_in_chunks(input_file):
    empty_line = '\n\n'
    with open(input_file, 'r') as reader:
        split_file = reader.read().split(empty_line)

    #    split_file = remove_empty_line(split_file)

    return remove_empty_strings(split_file)


def split_in_lines_and_remove_empty(definition):
    lines = definition.split('\n')
    for line in lines:
        if len(line) < 1:
            lines.remove(line)

    return lines


def convert_if_switch_parameter(name):
    if name[0] == '$':
        name = "switch_" + name[1:]
    return name


def parse_param_list(definition):
    # parameter definition starts in the second line
    param_lines = split_in_lines_and_remove_empty(definition)[1:]

    param_names = {}

    for line in param_lines:
        param_definition = line.split()

        if len(param_definition) >= 2:
            param = param_definition[0]

            # depending on if it is a switch parameter or not, the value is found in the 2nd or 3rd entry
            if len(param_definition) == 2:
                value = param_definition[1]
            else:
                value = param_definition[2]

            param = param.lower()
            param = convert_if_switch_parameter(param)
        else:
            print(" *** error ***: could not parse the parameter in the following model definition")
            print(definition)
            return None

        param_names[param] = value

    return param_names


#
# parse the lmodel.dat file to get all defined local model names
# input:  filename for the model definition (lmodel.dat file)
# return: list of model names
def get_model_names(lmodeldat_file):
    definition_list = []

    for single_definition_string in read_file_in_chunks(lmodeldat_file):
        definition_list.append(ModelDefinition(single_definition_string))

    return definition_list


def get_wrapper_lmod(local_model_name, function_name):
    parameter_list = "const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxError, const char *init"
    function_call = "xspec_C_wrapper_eval_model(ModelName::" + local_model_name + ", parameter, flux, Nflux, energy);"

    return f"""
extern "C" void {function_name}({parameter_list}) 
{{
    {function_call}
}} 
"""


def write_std_header(file):
    file.write(header)


def write_xspec_wrapper_file(outfile_name, model_definition):
    includes = "#include \"LocalModel.h\"\n"

    file = open(outfile_name, "w")

    write_std_header(file)
    file.write(includes)

    for model in model_definition:
        file.write(get_wrapper_lmod(model.model_name, model.function_name))

    file.close()


def get_implemented_lmod(local_model_name, param_list):
    param_class = "XPar::"

    lmodel_str = \
        f"""   {{ModelName::{local_model_name},
    XspecParamList(\"{local_model_name}\", 
        {{ """

    separator = ""
    for par_key in param_list.keys():
        lmodel_str += separator +param_class + par_key
        separator = ", "

    lmodel_str += "}, \n\t {"
    separator = ""
    for par_key in param_list.keys():
        lmodel_str += separator + param_list[par_key]
        separator = ", "


    lmodel_str += "})\n   },\n"

    return lmodel_str


def write_class_definition(file):
    class_definition_cpp = """

#include "ModelDefinition.h"

class XspecParamList: public ParamList {

 public:
  XspecParamList(std::string model_name, std::vector<XPar> parnames, std::vector<double> parvalues):
      ParamList(std::move(parnames), std::move(parvalues)) {
        m_name = std::move(model_name);
      };

  [[nodiscard]] std::string name() const {
    return m_name;
  }

 private:
  std::string m_name;
};


class XspecModelDatabase{

 public:
  std::string name_string(ModelName name) const{
    return model_definition.at(name).name();
  }

  ParamList param_list(ModelName name) const{
    return model_definition.at(name);
  }

  std::unordered_map<ModelName, XspecParamList> all_models() const{
    return model_definition;
  }

 private:
  const std::unordered_map<ModelName, XspecParamList> model_definition = {
  {"""
    file.write(class_definition_cpp)


def write_model_database(file, definition):
    for model in definition:
        file.write(get_implemented_lmod(model.model_name, model.params))
    file.write("  }\n };\n")


def write_xspec_implement_models(outfile_name, model_definition):
    includes = ""

    file = open(outfile_name, "w")

    write_std_header(file)

    c_define_string = "RELXILL_XSPEC_WRAPPER_H_"
    file.write(f"#ifndef {c_define_string}\n")
    file.write(f"#define {c_define_string}\n\n")

    file.write(includes)

    write_class_definition(file)

    write_model_database(file, model_definition)

    file.write("\n };\n")

    file.write(f"\n#endif //{c_define_string}")

    file.close()


if __name__ == '__main__':

    if len(sys.argv) == 3:
        input_lmodel_file = sys.argv[1]
        output_wrapper_file_cpp = sys.argv[2]
        output_wrapper_file_header = output_wrapper_file_cpp.replace(".cpp", ".h")

        print(
            f"\n *** creating {output_wrapper_file_cpp} and {output_wrapper_file_header} by parsing {input_lmodel_file}:")

        model_definition = get_model_names(input_lmodel_file)
        write_xspec_wrapper_file(output_wrapper_file_cpp, model_definition)
        write_xspec_implement_models(output_wrapper_file_header, model_definition)

    else:
        print("""
    Usage: ./create_lmod_wrapper.py [lmodel.dat] [xspec_wrapper.cpp]
        """)
        exit(1)
