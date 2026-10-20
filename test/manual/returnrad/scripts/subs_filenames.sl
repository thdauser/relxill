define clean_output_files(){
   () = system("rm -f __relxillOutput*");
}

define getFilenameStruct(){
   
   variable outdir = "./";
   variable fbase = outdir + "__relxillOutput-rrad-bbody-%s.fits";
   variable fbase_dat = outdir + "__relxillOutput-rrad-bbody-%s.dat";

   variable f_rframe = struct{
      specRet  = sprintf(fbase, "rframe-bbodyRet" ),      %% returning spectrum on the disk (rest-frame)
      specPri  = sprintf(fbase, "rframe-bbodyEmit" ),      %% primary BBody spectrum (rest-frame)
      xillPrim = sprintf(fbase, "rframe-primary" ),  %% emitted spectrum for boost=0 (rest-frame, identical to specPri)
      xillRefl = sprintf(fbase, "rframe-reflect" ),  %% reflected xillver spectrum, boost=-1 (rest-frame)
      xillCombined = sprintf(fbase, "rframe-total" )  %% reflected+primary (bbodyEmit) spectrum, boost=1 (rest-frame)
   }
   ;
   
   variable f_obs = struct{    
      reflect = sprintf(fbase, "obs-reflect" ),
      primary = sprintf(fbase, "obs-primary" ),
      total = sprintf(fbase, "obs-total" ),
      mirror  = sprintf(fbase, "obs-mirror-total" ),
      mirrorPrim  = sprintf(fbase, "obs-mirror-primary" ),
      mirrorRefl  = sprintf(fbase, "obs-mirror-reflect" )
   }
   ;
   
   
   variable f_disk = struct{
      rrad = sprintf(fbase_dat,"emis_profile")
   }
   ;
	 
   variable files = struct{
      rframe = f_rframe,
      fobs   = f_obs,
      disk   = f_disk
   }
   ;
   
   return files;
}


define read_ion_gradient(_fname){
%    rlo [R_g]      rhi [R_g]       log(xi)         log(N)       delta_emit F_x
   variable rlo, rhi, logxi, logn, del_emit, emis;
   (rlo, rhi, logxi, logn, del_emit, emis) = readcol(_fname,1,2,3,4,5,6);

   return struct{rlo=rlo, rhi=rhi, logxi=logxi, logn=logn, del_emit=del_emit, emis=emis};
}

define convert_plot_png(fname){
   () = system(sprintf("convert -density 200 %s %s",
   fname,
   strreplace(fname, ".pdf", ".png") )
   );
}

define render_with_png(pl, fname){
   pl.render(fname);
   () = system(sprintf("convert -density 200 %s %s",
   fname,
   strreplace(fname, ".pdf", ".png") )
   );
}


