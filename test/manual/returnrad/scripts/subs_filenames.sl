define clean_output_files(){
   () = system("rm -f __relxillOutput*");
}

define getFilenameStruct(){
   
   variable outdir = "./";
   variable fbase = outdir + "__relxillOutput-rrad-bbody-%s.fits";
   variable fbase_dat = outdir + "__relxillOutput-rrad-bbody-%s.dat";

   %% all files are fluxes; rest frame: angle-averaged flux per unit area of the zone (see the table in
   %% doc/returnrad_disk_emission.md, "Diagnostic output files")
   variable f_rframe = struct{
      specRet  = sprintf(fbase, "rframe-bbodyRet" ),      %% returning flux incident on the disk (rest-frame)
      specPri  = sprintf(fbase, "rframe-bbodyEmit" ),      %% primary BBody flux (rest-frame)
      xillRefl = sprintf(fbase, "rframe-reflect" ),  %% reflected xillver flux (angle-averaged), scaled with |boost| (rest-frame)
      mirrorRefl = sprintf(fbase, "rframe-mirror-reflect" )  %% mirrored returning flux, |boost| * specRet (rest-frame)
   }
   ;
   
   %% observed spectra of the zones (convolved with the line profile of each zone)
   variable f_obs = struct{    
      reflect = sprintf(fbase, "obs-reflect" ),          %% reflected xillver flux, scaled with |boost|
      primary = sprintf(fbase, "obs-bbodyEmit" ),        %% primary BBody flux (written for any boost)
      mirrorRefl  = sprintf(fbase, "obs-mirror-reflect" )  %% mirrored returning flux, scaled with |boost|
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


