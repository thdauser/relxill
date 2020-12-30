#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

require("isisscripts");
require("fits_model_struct");
require("subs_fitfunctions.sl");
require("test_setup.sl");


load_relxill_model_devel(modlib);

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

   __set_hard_limits("relxilllp","h",-100,1000);
   fit_fun("relxilllp");

   
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(filename_refdata));
   set_par("*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_refdata_filename(filename_refdata));
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_refdata_filename(filename_refdata));   
   set_par("*.fixReflFrac",2);
   fits_write_model_struct(get_refdata_filename(filename_refdata));      
   set_par("*.h",-1.1,0,-10,100);
   fits_write_model_struct(get_refdata_filename(filename_refdata));
   set_par("*.fixReflFrac",1);
   fits_write_model_struct(get_refdata_filename(filename_refdata));
}
%}}}

%%% MAIN %%%

() = system("rm -f $DATA_DIR/*/*refdat*"$);

variable filename_default = "%s_defparam_refdat_%04i.fits";
variable filename_random  = "%s_random_refdat_%04i.fits";
variable filename_special  = "%s_special_refdat_%04i.fits";
variable num_random_evaluations = 5;

variable ff;
foreach ff(get_implemented_fitfunctions()){
   create_default_refdata(ff, filename_default);
   create_random_refdata(ff, filename_random, num_random_evaluations);
}
create_refdata_relxilllp(filename_special);

%%%%%%%%%%