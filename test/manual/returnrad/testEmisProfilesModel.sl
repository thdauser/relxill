require("isisscripts");


load_xspec_local_models("build/");
putenv("RELXILL_WRITE_OUTFILES=1");


define readStdEmisProfile(fname){ %{{{
   
   
      variable rad, emis;
      (rad, emis) = readcol(fname,1,2);
   
      remove(fname);
   
      variable indAscending = array_sort(rad);
      rad  = rad[indAscending];
      emis = emis[indAscending];
   
      return struct{rad=rad, emis=emis};
}
%}}}


variable lo,hi;
(lo,hi) = log_grid(0.5,9.5,500);

%variable a = 0.998;
%variable rin = [kerr_rms(0.995), kerr_rms(0.991), kerr_rms(0.99), kerr_rms(0.989)];

variable a = [0.995,0.995,0.995,0.995];
variable rin = [kerr_rms(0.99025), kerr_rms(0.9902), kerr_rms(0.990), kerr_rms(0.9901)];
print(rin);

variable n = length(rin);

variable pl = xfig_plot_new();
pl.world(rin[0]*0.99,20,5e-4,2e-1;loglog);
variable plr = xfig_plot_new();
plr.world(rin[0]*0.99,20,0.95,1.05;xlog);

__set_hard_limits("rellinelpRet","Rin",-100,100);
fit_fun("rellinelpRet");

variable ii;
variable val = Array_Type[n];
variable emis_rr = Struct_Type[n];
variable emis_in = Struct_Type[n];

_for ii(0,n-1){
   set_par("*.a",a[ii]);
   set_par("*.Rin",rin[ii],0,1.2,100);
   set_par("*.return_rad",-1);
   set_par("*.h",3);
   set_par("*Rout",1000,0,100,1000);
   
   val[ii] = eval_fun_keV(lo,hi)/(hi-lo)*lo;
 
   emis_rr[ii] = readStdEmisProfile("test_emis_profile_rrad_output.dat");
   emis_in[ii] = readStdEmisProfile("test_emis_profile_rrad_input.dat");
}
   
variable rat_ind = 2;

_for ii(0,n-1){
   pl.plot(emis_rr[ii].rad,emis_rr[ii].emis;color=ii+1);
   pl.plot(emis_in[ii].rad,emis_in[ii].emis;color=ii+1);
   
   variable ipol_emis = interpol_points(emis_rr[rat_ind].rad,emis_rr[ii].rad,emis_rr[ii].emis);
   plr.plot(emis_rr[rat_ind].rad,ipol_emis/emis_rr[rat_ind].emis;color=ii+1,line=-1,sym="point");
}


xfig_multiplot(pl,plr).render("plots/testEmisProfilesIpol.pdf");