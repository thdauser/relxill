#!/usr/bin/env isis-script

_traceback=1;

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

require("create_output_files.sl");

variable ff_plots = "relxillBB";
variable prefix = ff_plots+"_";

%% =1= set parameters (currently the only ones for RTR table v0.3)
variable spin = 0.998; 
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
   variable plotTitle = sprintf("Rest frame, per zone normalized to $N_\mathrm{xill}$ of the return: Return (solid), Primary (dashed) -- spin=%.3f"R,
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.rframe.specRet,files.rframe.specPri];
   variable iRframe = [1,1];
   
   %% the primary is about 1/(returning fraction) ~ 25 times above the returning spectrum
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec, ymax=1e5);
   
   pl.render(dir_plots+prefix+"diagnosePlotRframe.pdf");
   
   return EXIT_SUCCESS;
}

define getDiagnoseXillverPlot(){
   
   variable plotTitle = sprintf("Rest Frame (per zone noramlized to return): Xillver (solid), Return (dashed) -- spin=%.3f"R,
				spin);

   message("\n"+plotTitle);

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
   
   variable plotTitle = sprintf("Observer, Mirror Refl: Reflected Mirror (solid), Primary BBody (dashed)  -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   %% (the order is needed for plot_kerrbb, which is normalized to the second file, the primary spectrum)
   variable fnames  = [files.fobs.mirrorRefl,files.fobs.primary];
   variable iRframe = [0,0];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec, plot_kerrbb);
   
   pl.render(dir_plots+prefix+"diagnosePlotMirrorBbody.pdf");
   
   return EXIT_SUCCESS;
}

define testDiagnosePlots(){

   variable params = get_param_struct(spin, Tin, logxi, logn);
   createOutputFiles(ff_plots,params);

   if (getDiagnoseMirrorBbody() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseRframePlot() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseXillverPlot() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseRelxillBBPlot() != EXIT_SUCCESS) return EXIT_FAILURE;

   return EXIT_SUCCESS;
}


if (length(__argv)>0){
   testDiagnosePlots();
}
