#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-


variable ff =     ["xillverCp", "xillver","xillverD","xillverNS","xillverCO"];
variable ff_tab = ["xillver-comp.fits", "xillver-a-Ec5.fits", "xillverD-5.fits", "xillverNS.fits","xillverCO.fits"];

variable ii, n = length(ff);

_for ii(0,n-1){
   message( sprintf("./create_xillspec_ref_single.sl %s %s",ff[ii],ff_tab[ii]) );
   () = system(sprintf("./create_xillspec_ref_single.sl %s %s",ff[ii],ff_tab[ii]) );
}