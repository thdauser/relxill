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
      xillPrim = sprintf(fbase, "rframe-xillverPrim" ),  %% primary xillver spectrum (rest-frame)
      xillRefl = sprintf(fbase, "rframe-xillverRefl" ),  %% reflected xillver spectrum (rest-frame)
      xillCombined = sprintf(fbase, "rframe-total" )  %% reflected+primary (bbodyEmit) spectrum (rest-frame)
   }
   ;
   
   variable f_obs = struct{    
      reflect = sprintf(fbase, "obs-reflect" ),
      primary = sprintf(fbase, "obs-primary" ),
      total = sprintf(fbase, "obs-total" ),
      mirror  = sprintf(fbase, "obs-mirror" ),
      mirrorPrim  = sprintf(fbase, "obs-mirror-primary" ),
      mirrorRefl  = sprintf(fbase, "obs-mirror-refl" )
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
