
define getFilenameStruct(){
   
   variable outdir = "build/";
   variable fbase = outdir + "debug-testrr-bbody-%s.fits";

   variable f_rframe = struct{
      specRet  = sprintf(fbase, "rframe-specRet" ),      %% returning spectrum on the disk (rest-frame)
      specPri  = sprintf(fbase, "rframe-specPri" ),      %% primary BBody spectrum (rest-frame)
      xillRefl = sprintf(fbase, "rframe-xillverRefl" ),  %% reflected xillver spectrum (rest-frame)
      xillPrim =  sprintf(fbase, "rframe-xillverPrim" )  %% primary spectrum producing reflected xillver spectrum (rest-frame)
   };
   %% Note that rframe_xillPrim is ideally the same as rframe_specRet,
   %% but here only trying to approximate this shape
   %% 
   variable f_obs = struct{    
      reflect = sprintf(fbase, "obs-reflect" ),
      primary = sprintf(fbase, "obs-primary" ),
      mirror  = sprintf(fbase, "obs-mirror" )
   };
   
   
   variable files = struct{
      rframe = f_rframe,
      fobs   = f_obs
   };
   
   return files;
}
