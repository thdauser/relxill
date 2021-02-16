#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

%%% call the routine like "./check_model_functions.sl DEV"
%%% if also the model in "lmodel_relxill_devel.dat" should 
%%% be tested
variable TEST_DEVEL = 0; 

if (__argc>1 && __argv[1]=="DEV"){
   TEST_DEVEL = 1;
}


_traceback=1;

%% we can load any model
variable ff = "relxilllpDCp";

variable modlib = "build/librelxill.so";

if (stat_file(modlib) == NULL){
   message("\n **** error : local relxill model not found ; exiting ... **** \n ");
   exit;
}

require("xspec");
load_xspec_local_models(modlib);
require("isisscripts");


putenv("DEBUG_RELXILL=1");
fit_fun("relxilllpDCp");
_traceback=1;


define fit_fun_default(ff){ %{{{
   %% only works for single model components %%
   
   fit_fun(ff);
   
   variable p,pa = get_params();   
   variable df;
   foreach p(pa){
      df = eval(sprintf("%s_default(%i)",qualifier("ff0",ff),p.index-1));
      p.value=df.value;
      p.min=df.min;
      p.max=df.max;
      p.freeze=df.freeze;
   }
   set_params(pa);
}
%}}}

variable lo0, hi0;
(lo0,hi0) = log_grid(0.5,1000,3000);


variable hbase = 3.0;
variable htop = [3.0,3.1,5.,8.,15.,20.,50,100.];


define read_emis_profile(){
   variable fname = "test_emis_profile.txt";
   variable rad, emis;
   (rad, emis) = readcol(fname,1,2);
   
   variable ind_sort = array_sort(rad);
   rad = rad[ind_sort];
   emis = emis[ind_sort];
   
   return struct{rad=rad, emis=emis};   
}



define read_relline_profile(){
   variable fname = "test_relline_profile.dat";
   variable ener, flux;
   (ener, flux) = readcol(fname,1,3);
      
   return struct{ener=ener, flux=flux};
}

define set_par_vals(hbase,htop,beta){

   fit_fun_default(ff);
   
   list_par;
   set_par("*.hbase",hbase);
   set_par("*.htop",htop);
   set_par("*.beta",beta);
   
   set_par("*.a",qualifier("a",0.998));

}


define eval_param(hbase,htop,beta){

   set_par_vals(hbase,htop,beta);

   () = eval_fun(1,2);
   
}

define eval_param_spec(hbase,htop,beta){

   set_par_vals(hbase,htop,beta);

   set_par("*.refl_frac",-1,0,-1,10);
   set_par("*.fixReflFrac",1);
   variable val = eval_fun_keV(lo0, hi0) / (hi0-lo0) * 0.5*(lo0+hi0);

   set_par("*.fixReflFrac",0);
   set_par("*.refl_frac",0,0,-1,10);
   variable val0 = eval_fun_keV(lo0, hi0) / (hi0-lo0) * 0.5*(lo0+hi0);

   return struct{bin_lo=lo0, bin_hi=hi0, flux=val, flux0=val0};
}


define get_emis_profile(hbase,htop,beta){
   
   eval_param(hbase,htop,beta);
   return read_emis_profile();
   
}

define get_relline_profile(hbase,htop,beta){
   
   eval_param(hbase,htop,beta);
   return read_relline_profile();
   
}


define get_std_plot(){
   
   variable pl = xfig_plot_new();
   
   pl.world(1.0,400,1e-7,1;loglog);
   pl.xlabel("radius");
   pl.ylabel("emissivity");
 
   return pl;
}


define get_std_plot_line(){
   
   variable pl = xfig_plot_new();
   
   pl.world(0.1,1,1e-7,1e-4;loglog);
   pl.xlabel("energy");
   pl.ylabel("flux [a.u.]");
 
   return pl;
}


define get_std_plot_spec(){
   
   variable pl = xfig_plot_new();
   
   pl.world(0.5,100,1e-1,1e3;loglog);
   pl.xlabel("energy");
   pl.ylabel("flux [a.u.]");
 
   return pl;
}



define plot_scale_htop(hbase,htop){
   
   variable beta = 0.0;
   
   variable ii, n = length(htop);
   variable emis = Struct_Type[n];
   variable relflux = Struct_Type[n];
   variable spec = Struct_Type[n];
   
   _for ii(0,n-1){
      emis[ii] = get_emis_profile(hbase,htop[ii],beta);
      relflux[ii] = get_relline_profile(hbase,htop[ii],beta);
      spec[ii] = eval_param_spec(hbase,htop[ii],beta);
   }
   
   
   variable pl = get_std_plot();
   variable plr = get_std_plot_line();
   variable pls = get_std_plot_spec();
   
   _for ii(0,n-1){
      pl.plot(emis[ii].rad, emis[ii].emis; color=CB_COLOR_SCHEME_NB[ii], width=3);

      plr.plot(relflux[ii].ener, relflux[ii].flux; color=CB_COLOR_SCHEME_NB[ii], width=3);

      pls.plot(spec[ii].bin_lo, spec[ii].flux; color=CB_COLOR_SCHEME_NB[ii], width=3);      
      pls.plot(spec[ii].bin_lo, spec[ii].flux0; color=CB_COLOR_SCHEME_NB[ii], width=3, line=1);
   }
   
   
   variable plsum = xfig_new_vbox_compound(pl,plr);
   xfig_new_vbox_compound(plsum, pls).render("plot_scale_htop.pdf");

}



plot_scale_htop(hbase, htop);