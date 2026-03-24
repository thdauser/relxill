#!/usr/bin/env isis-script

_traceback=1;

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

require("create_output_files.sl");

variable ff_plots = "relxillBBxill";
variable prefix = ff_plots+"_";

%% =1= set parameters
variable spin = 0.998;  %% currently only done for this value
variable Tin = 2.0;
variable logxi = 2.0;
variable logn = 18.0;



define eval_std_model(){
%{{{

   variable lo,hi;
   (lo, hi) = log_grid(0.1,40,3000);

   variable mn = 0.5*(lo+hi);

   set_par("*.boost",-1,0,-1,10);
   variable refl = eval_fun_keV(lo,hi)/(hi-lo)*mn^2;

   set_par("*.boost",0);
   variable direct = eval_fun_keV(lo,hi)/(hi-lo)*mn^2;

   set_par("*.boost",1);
   variable full = eval_fun_keV(lo,hi)/(hi-lo)*mn^2;


   return struct{lo=lo,mn=mn,hi=hi, refl=refl, direct=direct, full=full};
}
%}}}

define getDiagnoseRframePlot(){
   variable plotTitle = sprintf("Restframe Plot of Return (solid) and Primary (dashed) -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.rframe.specRet,files.rframe.specPri];
   variable iRframe = [1,1];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle);
   
   pl.render(dir_plots+prefix+"diagnosePlotRframe.pdf");
   
   return EXIT_SUCCESS;
}

define getDiagnoseXillverPrimPlot(){
   
   variable plotTitle = sprintf("Plot of Real Primary  (solid) and Xillver Primary (dashed) used for the reflection -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.rframe.specRet,files.rframe.specPri];
   variable iRframe = [1,1];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec);
   
   pl.render(dir_plots+prefix+"diagnosePlotXillverPrim.pdf");
   
   return EXIT_SUCCESS;
}

define getDiagnoseXillverPlot(){
   
   variable plotTitle = sprintf("Rest-frame: Xillver (solid), and Xillver Prim (dashed) -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.rframe.xillRefl,files.rframe.specRet];
   variable iRframe = [1,1];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec);
   
   pl.render(dir_plots+prefix+"diagnosePlotXillver.pdf");
   
   return EXIT_SUCCESS;
}

define getDiagnoseRelxillBBPlot(){
%{{{
   
   variable plotTitle = sprintf("Relat Refl (solid), Primary BBody (dashed), Mirrord BBody (dotted)",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.fobs.reflect,files.fobs.primary, files.fobs.mirrorRefl];
   variable iRframe = [0,0,0];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle);
   
   pl.render(dir_plots+prefix+"diagnosePlotRelxillBB.pdf");
   
   return EXIT_SUCCESS;
}
%}}}

define getDiagnoseMirrorBbody(){
%{{{
   
   variable plotTitle = sprintf("Observer, Mirror Refl: Combined (solid), Primary BBody (dashed), Reflected Mirror (dotted)  -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.fobs.mirror,files.fobs.mirrorPrim,files.fobs.mirrorRefl];
   variable iRframe = [0,0,0];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec, plot_kerrbb);
   
   pl.render(dir_plots+prefix+"diagnosePlotMirrorBbody.pdf");
   
   return EXIT_SUCCESS;
}
%}}}


define getCompareXillverTables(){

   variable params = get_param_struct(spin, Tin, logxi, logn);
   
   load_std_model("relxillBBxill",params);
   variable rtr = eval_std_model();
   
   load_std_model("relxillBB",params);
   variable xillns = eval_std_model();

   variable pl = xfig_plot_new();
   pl.world(0.1,40,1,900;loglog);

   pl.plot(rtr.mn,rtr.refl;width=2,line=0,color="red");
   pl.plot(rtr.mn,rtr.direct;width=2,line=1);
   pl.plot(rtr.mn,rtr.full;width=2,line=1,color="red");

   pl.plot(xillns.mn,xillns.refl;width=2,line=0,color="blue");
   pl.plot(xillns.mn,xillns.direct;width=2,line=1);
   pl.plot(xillns.mn,xillns.full;width=2,line=1,color="blue");

   
   pl.render(dir_plots+"diagnosePlotTables.pdf");

   return EXIT_SUCCESS;
}

define testDiagnosePlots(){

   variable params = get_param_struct(spin, Tin, logxi, logn);
   createOutputFiles(ff_plots,params);

   if (getDiagnoseMirrorBbody() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseRframePlot() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseXillverPlot() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseRelxillBBPlot() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getCompareXillverTables() != EXIT_SUCCESS) return EXIT_FAILURE;

   return EXIT_SUCCESS;
}


if (length(__argv)>0){
   testDiagnosePlots();
}
