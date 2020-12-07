#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-


require("isisscripts");
require("fits_model_struct");
require("subs_fitfunctions.sl");
require("test_setup.sl");

variable DATA_DIR = "refdata_localModels/";

() = system("rm -f $DATA_DIR/refdata*"$);


load_relxill_model_devel(modlib);

_traceback=1;

variable counter = 0;
variable std_fname = "refdata_%04i_%s_mod.fits";
variable PAR_DIR = "parfiles/";

define get_fname(fname){ %{{{
   counter++;
   variable funname = qualifier("par",get_fit_fun());
   return DATA_DIR+sprintf(fname,counter,funname);
}
%}}}

define test_xillver(){ %{{{
   fit_fun("xillver");
   set_par("*.refl_frac",-1,0,-10,10);
   set_par("*.Incl",30.1);
   fits_write_model_struct(get_fname(std_fname));

   set_par("*.logxi",0);
   fits_write_model_struct(get_fname(std_fname));

   set_par("*.logxi",3.1);
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));

   set_par("*.refl_frac",1,0,-10,10);

   variable elo, ehi;
   (elo,ehi) = log_grid(3,200,500);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,50,500);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,200,2000);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   
}
%}}}
define test_relxillxillver(){ %{{{
   fit_fun("relxill+xillver");
   set_par("xillver*.refl_frac",-1);
   set_par("relxill*.refl_frac",0.3,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));

   set_par("relxill*.refl_frac",1,0,-10,10);
   variable elo, ehi;
   (elo,ehi) = log_grid(3,200,500);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,50,500);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,200,2000);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   set_par("relxill*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));


}
%}}}
define test_relxilllpxillver(){ %{{{
   fit_fun("relxilllp+xillver");
   set_par("relxilllp*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("relxilllp*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("relxilllp*.refl_frac",0);
   fits_write_model_struct(get_fname(std_fname));
   set_par("relxilllp*.refl_frac",1);
   set_par("relxilllp*.fixReflFrac",2);
   fits_write_model_struct(get_fname(std_fname));      
   set_par("xillver*.refl_frac",-0.3,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));   
}
%}}}
define test_relxill(){ %{{{
   fit_fun("relxill");
   set_par("*.refl_frac",1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",0);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.refl_frac",-1,0,-10,10);
   fits_write_model_struct(get_fname(std_fname));
   
   set_par("*.logxi",0);
   fits_write_model_struct(get_fname(std_fname));

   set_par("*.logxi",3.1);
   set_par("*.refl_frac",-1,0,-10,10);
   set_par("*.z",0.1);
   fits_write_model_struct(get_fname(std_fname));

   set_par("*.z",0.5);
   fits_write_model_struct(get_fname(std_fname));
   set_par("*.z",0.0);

   
   set_par("*.refl_frac",1,0,-10,10);


   variable elo, ehi;
   (elo,ehi) = log_grid(3,200,500);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,50,500);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);

   (elo,ehi) = log_grid(3,200,2000);
   fits_write_model_struct(get_fname(std_fname);elo=elo,ehi=ehi);


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


define test_parfiles(){ %{{{
   
   variable pars = glob(PAR_DIR+"*.par");
   
   variable par;
   variable parname;
   foreach par(pars){
      load_par(par);
      parname = strchop(par)[0];
      fits_write_model_struct(get_fname(std_fname;par=parname));      
   }
   
}
%}}}



define create_defaultRefdata(ff, filename_refdata){ %{{{
 
   fit_fun(ff);

   variable outfile = get_fname(filename_refdata);
   vmessage(" - creating refdata: %s ", outfile);
   fits_write_model_struct(outfile);
   
}
%}}}
define create_defaultRefdata_allFitFunctions(){ %{{{
   
   variable filename_refdata = "refdat_defparam_%04i_%s.fits";
   
   variable ff;
   foreach ff(get_implemented_fitfunctions()){
      create_defaultRefdata(ff, filename_refdata);
   }
}
%}}}


create_defaultRefdata_allFitFunctions();



%test_xillver();
%test_relxill();
%test_relxillxillver();

%test_relxilllp;
