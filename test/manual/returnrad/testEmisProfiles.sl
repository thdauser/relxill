require("isisscripts");

variable outdir = "build/";
variable files = outdir+["test_emis_profile.dat","test_emis_profile_returnrad.dat","test_emis_profile_input.dat"];
variable filesIpol = outdir+["test_emisCoarse.dat","test_emisFineReference.dat","test_emisRebin.dat"];
variable col = CB_COLOR_SCHEME_NB;

define read_emis_profile(fname){
   variable rad, emis;
   (rad, emis) = readcol(fname, 1,2);
   return struct{rad=rad, emis=emis};
}

define plot_single_zones(pl){
   variable files_single_zone = outdir + "test_rrad_single_zone_"+["00","12","24","48"]+".dat";
     
   variable ii, n = length(files_single_zone);
   
   _for ii(0,n-1){
      variable dat = read_emis_profile(files_single_zone[ii]);
      pl.plot(dat.rad,dat.emis;line=ii);
   }
     
   
}


define plotEmisProfiles(files){

   variable pl = xfig_plot_new();
   
   variable ii, n = length(files);
   variable emis = Struct_Type[n];
   
   _for ii(0,n-1){
      
      emis[ii] = read_emis_profile(files[ii]);
   
      pl.world(1,1000,1e-12,100;loglog);
      pl.plot(emis[ii].rad,emis[ii].emis; width=2, line=ii, color=col[ii], depth=-ii);
   
   }
   
   plot_single_zones(pl);
   
   return pl;
}




plotEmisProfiles(files).render("plots/testEmisProfiles.pdf");

plotEmisProfiles(filesIpol).render("plots/testEmisProfilesIpol.pdf");