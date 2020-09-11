require("isisscripts");

variable outdir = "build/";
variable files = outdir+["test_emisRebin.dat","test_emisFineReference.dat"];
variable col = CB_COLOR_SCHEME_NB;

define read_emis_profile(fname){
   variable rad, emis;
   (rad, emis) = readcol(fname, 1,2);
   return struct{rad=rad, emis=emis};
}

variable pl = xfig_plot_new();

variable ii, n = length(files);
variable emis = Struct_Type[n];

_for ii(0,n-1){
   
   emis[ii] = read_emis_profile(files[ii]);
   
   pl.world(1,1000,1e-13,0.1;loglog);
   pl.plot(emis[ii].rad,emis[ii].emis; width=2, color=col[ii], line=ii, depth=-ii);
   
}


pl.render("plots/testEmisProfiles.pdf");