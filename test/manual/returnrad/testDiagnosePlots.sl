% -*- mode: slang; mode: fold -*-

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");


variable spin = 0.998;  %% currently only done for this value
variable Tin = 0.7;

define getDiagnoseRframePlot(){ %{{{
   
   variable plotTitle = sprintf("Restframe Plot of Return (solid) and Primary (dashed) -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.rframe.specRet,files.rframe.specPri];
   variable iRframe = [1,1];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle);
   
   pl.render(dir_plots+"diagnosePlotRframe.pdf");
   
   return EXIT_SUCCESS;
}
%}}}


define getDiagnoseXillverPrimPlot(){ %{{{
   
   variable plotTitle = sprintf("Plot of Real Primary  (solid) and Xillver Primary (dashed) used for the reflection -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.rframe.specRet,files.rframe.xillPrim];
   variable iRframe = [1,1];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec);
   
   pl.render(dir_plots+"diagnosePlotXillverPrim.pdf");
   
   return EXIT_SUCCESS;
}
%}}}


define getDiagnoseXillverPlot(){ %{{{
   
   variable plotTitle = sprintf("Rest-frame: Xillver (solid), and Xillver Prim (dashed) -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.rframe.xillRefl,files.rframe.xillPrim,files.rframe.specRet];
   variable iRframe = [1,1,1];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec);
   
   pl.render(dir_plots+"diagnosePlotXillver.pdf");
   
   return EXIT_SUCCESS;
}
%}}}





define getDiagnoseRelxillBBPlot(){ %{{{
   
   variable plotTitle = sprintf("Relat Refl (solid), Primary BBody (dashed), Mirrord BBody (dotted)",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.fobs.reflect,files.fobs.primary, files.fobs.mirrorRefl];
   variable iRframe = [0,0,0];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle);
   
   pl.render(dir_plots+"diagnosePlotRelxillBB.pdf");
   
   return EXIT_SUCCESS;
}
%}}}



define getDiagnoseMirrorBbody(){ %{{{
   
   variable plotTitle = sprintf("Observer, Mirror Refl: Combined (solid), Primary BBody (dashed), Reflected Mirror (dotted)  -- spin=%.3f",
				spin);


   variable files = getFilenameStruct();
   variable fnames  = [files.fobs.mirror,files.fobs.mirrorPrim,files.fobs.mirrorRefl];
   variable iRframe = [0,0,0];
   
   variable pl = get_2d_plot(fnames, iRframe, spin; title=plotTitle, noSumSpec, plot_kerrbb);
   
   pl.render(dir_plots+"diagnosePlotMirrorBbody.pdf");
   
   return EXIT_SUCCESS;
}
%}}}



define createOutputFiles(){ %{{{
   load_xspec_local_models("build/");
   fit_fun("relxillBB");
   
   putenv("RELXILL_WRITE_FILES=1");

   () = system("rm -f debug-testrr-bbody*");
   
   message("  Creating Output Files (calling relxillBB model) ");
   
   set_par("*.a",spin);
   set_par("*.kTbb", Tin);
%   () = eval_fun(1,2);
   
   %% creating the "mirror" files
   putenv("RELXILL_BBRET_NOREFL=1");

   set_par("*.Rout",400);
   set_par("*.boost",1.0);
   set_par("*.Incl",40);
   set_par("*.logxi",1.0);
   set_par("*.boost",1);
   () = eval_fun(1,2);

   set_par("*.boost",0);
   () = eval_fun(1,2);

   set_par("*.boost",-1.0,0,-10,10);
   () = eval_fun(1,2);
   

   %% creating the "reflection" files
   putenv("RELXILL_BBRET_NOREFL=0");

   set_par("*.boost",1);
   () = eval_fun(1,2);
   
   set_par("*.boost",-1.0,0,-10,10);
   () = eval_fun(1,2);

   set_par("*.boost",0.0,0,-10,10);
   () = eval_fun(1,2);

   putenv("RELXILL_WRITE_FILES=0");
   putenv("RELXILL_BBRET_NOREFL=0");
}
%}}}



define testDiagnosePlots(){

   createOutputFiles();

   if (getDiagnoseMirrorBbody() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseRframePlot() != EXIT_SUCCESS) return EXIT_FAILURE;
   if (getDiagnoseXillverPlot() != EXIT_SUCCESS) return EXIT_FAILURE;
   
   if (getDiagnoseRelxillBBPlot() != EXIT_SUCCESS) return EXIT_FAILURE;

   return EXIT_SUCCESS;
}


if (length(__argv)>0){
   testDiagnosePlots();
}
