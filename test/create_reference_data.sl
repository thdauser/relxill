#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-


%% we can load any model
if (length(__argv)>1){
   load_xspec_local_models(__argv[1]);
}

require("isisscripts");
require("fits_model_struct");

variable counter = 0;
variable std_fname = "test_%s_%04i.fits";
variable DATA_DIR = "refdata/";

define get_fname(fname){ %{{{
   counter++;
   return DATA_DIR+sprintf(fname,get_fit_fun(),counter);
}
%}}}

define test_relxill(){ %{{{
   fit_fun("relxill");
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_fname(std_fname));
}
%}}}
define test_xillver(){ %{{{
   fit_fun("xillver");
   set_par("*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
}
%}}}
define test_relxillxillver(){ %{{{
   fit_fun("relxill+xillver");
   set_par("relxill*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("relxill*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("relxill*.refl_frac",0);
   fits_write_model_struct(get_fname(std_fname));
}
%}}}
define test_relxilllp(){ %{{{

   __set_hard_limits("relxilllp","h",-100,1000);
   fit_fun("relxilllp");

   
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_fname(std_fname));   
   set_par("*.fixReflFrac",2);
   fits_write_model_struct(get_fname(std_fname));      
   set_par("*.h",-1.1,0,-10,100);
   fits_write_model_struct(get_fname(std_fname));
}
%}}}

test_relxilllp;
