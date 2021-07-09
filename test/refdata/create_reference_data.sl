#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

require("isisscripts");
require("fits_model_struct");
require("subs_fitfunctions.sl");
require("test_setup.sl");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Crete reference data (from the model built in build/) %%%
%
% Usage: %  ./test_refdata_relxill.sl  [<model_name>]
% 
% If called without arguments, ALL fit functions are called. If 
% a model name is given as argument, only reference data for this 
% model is created.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variable dry_run = 0;

% putenv("RELXILL_WRITE_OUTFILES=1");

load_relxill_model_devel(modlib);

define get_fit_functions(){ %{{{
   if (length(__argv)>1){ 
      return __argv[1]; 
   } else {
      return get_implemented_fitfunctions(;; __qualifiers() );
   }
}
%}}}



variable DATA_DIR = "refdata_localModels/";

_traceback=1;
variable counter = 0;
variable std_fname = "refdata_%04i_%s_mod.fits";
variable PAR_DIR = "parfiles/";

define get_refdata_filename(fname){ %{{{
   counter++;
   if (counter>9999){
      message(" *** error: counter overflow ");
       exit(1);
   }
   
   variable local_model_name = qualifier("par",get_fit_fun());
   variable subdir = DATA_DIR + local_model_name + "/";
   
   () = system("mkdir -p $subdir"$);
   
   return subdir+sprintf(fname,local_model_name,counter);
}
%}}}

define test_xillver(){ %{{{
   fit_fun("xillver");
   set_par("*.refl_frac",-1,0,-10,10);
   set_par("*.Incl",30.1);
   fits_write_model_struct(get_refdata_filename(std_fname));

   set_par("*.logxi",0);
   fits_write_model_struct(get_refdata_filename(std_fname));

   set_par("*.logxi",3.1);
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_refdata_filename(std_fname));
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));

   set_par("*.refl_frac",1,0,-10,10);

   variable elo, ehi;
   (elo,ehi) = log_grid(3,200,500);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,50,500);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,200,2000);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   
}
%}}}
define test_relxillxillver(){ %{{{
   fit_fun("relxill+xillver");
   set_par("xillver*.refl_frac",-1);
   set_par("relxill*.refl_frac",0.3,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));

   set_par("relxill*.refl_frac",1,0,-10,10);
   variable elo, ehi;
   (elo,ehi) = log_grid(3,200,500);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,50,500);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,200,2000);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   set_par("relxill*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));


}
%}}}
define test_relxilllpxillver(){ %{{{
   fit_fun("relxilllp+xillver");
   set_par("relxilllp*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));
   set_par("relxilllp*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));
   set_par("relxilllp*.refl_frac",0);
   fits_write_model_struct(get_refdata_filename(std_fname));
   set_par("relxilllp*.refl_frac",1);
   set_par("relxilllp*.fixReflFrac",2);
   fits_write_model_struct(get_refdata_filename(std_fname));      
   set_par("xillver*.refl_frac",-0.3,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));   
}
%}}}
define test_relxill(){ %{{{
   fit_fun("relxill");
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_refdata_filename(std_fname));
   set_par("*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(std_fname));
   
   set_par("*.logxi",0);
   fits_write_model_struct(get_refdata_filename(std_fname));

   set_par("*.logxi",3.1);
   set_par("*.refl_frac",-1,0,-10,10);
   set_par("*.z",0.1);
   fits_write_model_struct(get_refdata_filename(std_fname));

   set_par("*.z",0.5);
   fits_write_model_struct(get_refdata_filename(std_fname));
   set_par("*.z",0.0);

   
   set_par("*.refl_frac",1,0,-10,10);


   variable elo, ehi;
   (elo,ehi) = log_grid(3,200,500);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,50,500);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,200,2000);
   fits_write_model_struct(get_refdata_filename(std_fname);elo=elo,ehi=ehi);


}
%}}}
define test_parfiles(){ %{{{
   
   variable pars = glob(PAR_DIR+"*.par");
   
   variable par;
   variable parname;
   foreach par(pars){
      load_par(par);
      parname = strchop(par)[0];
      fits_write_model_struct(get_refdata_filename(std_fname;par=parname));      
   }
   
}
%}}}
private define write_all_refl_frac_combinations(filename_refdata){ %{{{

   set_par("*.boost",1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(filename_refdata);; __qualifiers() );
   set_par("*.boost",0,0,-10,10);
   fits_write_model_struct(get_refdata_filename(filename_refdata);; __qualifiers() );
   set_par("*.boost",-1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(filename_refdata);; __qualifiers() );   
   set_par("*.boost",1);
   
}
%}}}


define create_default_refdata(ff, filename_refdata){ %{{{
 
   fit_fun(ff);
   fits_write_model_struct(get_refdata_filename(filename_refdata); verbose);
}
%}}}
define create_random_refdata(ff, filename_refdata, num_random){ %{{{
 
   fit_fun(ff);

   thaw("*.Incl");  
   freeze("*.norm"); %% we do not want to simply change the norm

   variable ii;
   _for ii(0, num_random-1){
      randomize_params();
      set_par("*.z",rand_flat(0.0,0.1)); %% special case, as typically frozen
      fits_write_model_struct(get_refdata_filename(filename_refdata); verbose);      
   }
}
%}}}
define create_refdata_relxilllp(filename_refdata){ %{{{

   message("\n creating special relxilllp reference data \n");
   __set_hard_limits("relxilllp","h",-100,1000);
   fit_fun("relxilllp");
   
   set_par("*.a",0.998);

   set_par("*.h",2.0,0,-100,100);
   write_all_refl_frac_combinations(filename_refdata; verbose);

   set_par("*.h",30,0,-100,100);
   write_all_refl_frac_combinations(filename_refdata; verbose);

   set_par("*.h",-1.1,0,-10,100);
   write_all_refl_frac_combinations(filename_refdata; verbose);
}
%}}}

%%% MAIN %%%

variable filename_default = "%s_defparam_refdat_%04i.fits";
variable filename_random  = "%s_random_refdat_%04i.fits";
variable filename_special  = "%s_special_refdat_%04i.fits";
variable num_random_evaluations = 3;

variable all_fit_functions = get_fit_functions();

variable ff;
foreach ff([all_fit_functions]){
   counter = 0;
   message(ff);
   if (dry_run==0){
      create_default_refdata(ff, filename_default);
      create_random_refdata(ff, filename_random, num_random_evaluations);
   }
}

if ( length(where(all_fit_functions=="relxilllp") ) > 0 ){
   message("relxillp - special");
   if (dry_run==0){
      counter = 0;
      create_refdata_relxilllp(filename_special);
   }
}

%%%%%%%%%%