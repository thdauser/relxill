require("isisscripts");
%require("scripts/subs_testSetup.sl");
%require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

variable spin = 0.998;
variable logxi = 2.0;

define model_eval(){

   putenv("RELXILL_WRITE_FILES=1");
   load_xspec_local_models("./build/");
   variable ff = "relxillBBxill";
   variable spin = 0.998;

   fit_fun(ff);
   set_par("*.a",spin);
   set_par("*.logxi",logxi);

   () = eval_fun(1,2);

   return spin;
}

%define get_fnames(){
%   variable f = getFilenameStruct();
%   return struct{prim=f.rframe.specPri , refl=f.rframe.specRet, emis=f.f_disk.rrad};
%}

variable rmax = 40;

define get_plot_iongrad(iongrad){

   variable pl = xfig_plot_new();

   pl.world(1,rmax,0,4;xlog);
   pl.world2(1,rmax,15,20;xlog);
   
   variable rmean = 0.5*(iongrad.rlo+iongrad.rhi);
   pl.plot(rmean, iongrad.logxi);
   pl.plot(rmean, iongrad.logn;world2,color="red");
   
   pl.ylabel("$\log \xi$"R);
   
   pl.y2axis(;color="red");
   pl.y2label("\color{red} $\log N$"R);
   
   pl.x2axis(;ticlabels=0);
   
   pl.xlabel("Radius $[R_\mathrm{g}]$"R);
   
   pl.title(sprintf(" Ion Gradient (a=%.3f - logxi=%.1f) ", spin, logxi));

   return pl;
}


define get_plot_emis_rrad(iongrad){
   variable pl = xfig_plot_new();

   pl.world(1,rmax,1e-9,1e-4;loglog);
   
   variable rmean = 0.5*(iongrad.rlo+iongrad.rhi);

   %% canonocal emissivity, normalized at the last bin
   variable emis_can = rmean^-3;
   variable ind = wherefirst(rmean >= 10.0);
   emis_can *= iongrad.emis[ind]/ emis_can[ind];
   

   pl.plot(rmean, iongrad.emis;color="blue");
   pl.plot(rmean, emis_can;color="gray");

   pl.ylabel("Emissivity");
   pl.xlabel("Radius $[R_\mathrm{g}]$"R);
   

   return pl;
}


model_eval();

variable f = getFilenameStruct();
variable iongrad = read_ion_gradient(f.disk.rrad);



variable pl = get_plot_iongrad(iongrad);

variable pl_emis = get_plot_emis_rrad(iongrad);

xfig_multiplot(pl, pl_emis).render("plots/iongradient_rrad.pdf");

