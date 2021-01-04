#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test the current model in build/ against the reference data %%%
%
% Usage: %  ./test_refdata_relxill.sl  [<model_name> [defparam|random] [update] ]
% 
% Note, that <model_name> can be a regular expression. For example, to
% the defparm test for all models, it would be:
%
%  ./test_refdata_relxill.sl '*' defparam
% 
% If 'update' is given as a last argument, the refdata will get updated 
% after the test. 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


require("isisscripts");
require("subs_fitfunctions.sl");
require("fits_model_struct");

require("test_setup.sl");

load_relxill_model_devel("build/librelxill.so");

define get_refdata_files(){ %{{{
   %% we can load any model
   variable ff = "*";
   if (length(__argv)>1){ ff = __argv[1]; }
   
   variable refdata_type = "*";
   if (length(__argv)>2){
      refdata_type = __argv[2]; 
  }
   
   variable glob_str = sprintf("%s/*_%s_ref*_*.fits",ff,refdata_type);
   variable fnames = glob(LMOD_REFDATA_DIR+sprintf(glob_str));
   
   return fnames[array_sort(fnames)];   
}
%}}}

define should_refdata_get_updated(){ %{{{
   if (__argv[-1] == "update"){
      return 1;
   } else {
      return 0;
   }
}
%}}}

define conv2tex(str){ %{{{
   return strreplace(str,`_`,`\_`);
}
%}}}

define add_labels(pl,goodness,fn){ %{{{

   variable ff_tex = conv2tex(get_fit_fun());
   variable str_tex = conv2tex(fn); 
   pl.add_object(xfig_new_text(sprintf(`%s [goodness %.3e] --- %s`,
				       ff_tex,goodness,str_tex	       )),
		 0.5,1.05,0,0;world0);
   
   variable params = get_params();
   
   variable ii, n=length(params);
   
   variable str;
   _for ii(0,n-1,1){
      str = sprintf(`%s : %.3e`, strreplace(params[ii].name,`_`,`\_`),
		    params[ii].value);
      pl.add_object(xfig_new_text(str;
				 size=`\footnotesize`),
		    1.02,0.95,-0.5,0+0.8*ii;world0);            
   }
}
%}}}
define plot_model_comparison(fn,dat,goodness){ %{{{
   
   variable pl = xfig_plot_new(15,8);
   variable plr = xfig_plot_new(15,4);
   
   variable elo = dat.bin_lo;
   variable ehi = dat.bin_hi;
   variable m_dat = dat.model;
   
   pl.world(min(elo),max(elo),0.1,200;loglog);
   plr.world(min(elo),max(elo),0.989,1.011;loglog);
   
   
   pl.hplot([elo,ehi[-1]], dat.value/(ehi-elo)*(0.5*(elo+ehi))^2;
	    loglog);
   pl.hplot([elo,ehi[-1]], m_dat/(ehi-elo)*(0.5*(elo+ehi))^2 ;
	    loglog, color="red", depth=-1 );
   
   plr.hplot([elo,ehi[-1]], dat.value*0+1; loglog);
   
   plr.hplot([elo,ehi[-1]], m_dat/dat.value; loglog, color="red");
   
   pl.ylabel(" Energy Flux ");
   plr.ylabel(" Ratio New/Old ");
   plr.xlabel("Energy [keV]" );
   
        
   add_labels(pl,goodness,fn);
   
   
   xfig_multiplot(pl,plr).render(fn+".pdf");
}
%}}}
private define fits_write_model_data(fn, dat){ %{{{
   variable fname_post = "_comparison.tmp";
   variable filename = fn+fname_post;
   fits_write_binary_table(filename, "MODEL", dat);
   vmessage(" evaluted model stored at %s ", filename);
}
%}}}



define check_single_model(fn){ %{{{

   variable status = EXIT_SUCCESS;

   variable dat = fits_read_model_struct(fn);
   variable m_dat = eval_fun_keV(dat.bin_lo,dat.bin_hi);
   dat = struct_combine(dat, struct{model=m_dat} );
   
   struct_filter(dat, where(dat.value != 0) );
   
   variable ind = where(dat.bin_lo>0.2);
   variable goodness =  sqrt(sum((dat.model[ind]/dat.value[ind]-1)^2))/length(dat.model[ind]);
   
   () = printf(" %s  \t deviation:  %.3e ",fn, goodness);   %% GOODNESS = sqr-distance / num_bins 
   
   if (goodness>goodness_limit_modelcomparison){
      status = EXIT_FAILURE;            
      message("    *** test FAILED *** ");
   } else {
      message("");
   }
   
   
   if ( status!=EXIT_SUCCESS  || qualifier_exists("plot") ){
      
      save_par(fn+".par");
      
      plot_model_comparison(fn,dat,goodness);
      fits_write_model_data(fn, dat);
   }
   
   if ( qualifier("update") == 1 ){
      vmessage(" *** UPDATING refdata: %s", fn);
      fits_write_model_struct(fn);
   }

      
   return status;
}
%}}}

define print_summary(status){ %{{{
   
   variable ntests = length(status);
   variable nfailed = length( where(status!=EXIT_SUCCESS) );
   
   message( " ******************************* ");
   vmessage(  " %i / %i tests successful ", ntests-nfailed, ntests);
   if (nfailed==0){      
      vmessage(  "    => SUCCESS ");
   } else {
      message(   "    => FAILED "); 
   }
   message( " ******************************* ");
   
}
%}}}



%%%% MAIN %%%% 

variable fnames = get_refdata_files();
variable update_refdata = should_refdata_get_updated();

variable ii, n = length(fnames);
variable status = Int_Type[n];

_for ii(0,n-1,1){
   status[ii] = check_single_model(fnames[ii]; update=update_refdata );   
}

print_summary(status);


if (length( where(status!=EXIT_SUCCESS)) > 0 ){
   exit(1);
}


%%%%%%%%%%%%%%
