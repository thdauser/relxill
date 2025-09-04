require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

define runtest(ffs){ 
      
   variable status = EXIT_SUCCESS; 
   
   variable fit_fun_rrad = "relxillBB";

   if (0 == is_defined(fit_fun_rrad) )     {
      msg_log += sprintf( "skipping test as function %s does not exist",fit_fun_rrad);
      vmessage( "     -> skipping test as function %s does not exist",fit_fun_rrad);
      return status;
   }
   
	   
   fit_fun(fit_fun_rrad);

   putenv("RELXILL_BBRET_NOREFL=0");

   variable Tin = 0.5;

   set_par("*.kTbb",Tin);
   set_par("*boost",0);
   set_par("*Rin",-10);
%   set_par("*Afe",5);
%   set_par("*a",0.95);
   set_par("*Incl",40);
   
%   list_par;
   
   variable lo, hi;
   (lo,hi) = log_grid(0.1,20,1000);
   variable val = eval_fun_keV(lo,hi) / (hi-lo)*lo;

   fit_fun("diskbb");
   set_par("*Tin",Tin);
   variable val_ref = eval_fun_keV(lo,hi) / (hi-lo)*lo;
   val_ref *= 1./val_ref[500]*val[500];

   variable ind = where(Tin < lo < 5*Tin);
   variable dE = (hi-lo);
   variable diff =  abs(1 - sum(val[ind]*dE[ind]) / sum(val_ref[ind]*dE[ind]));

   %% profiles will not be the same as
   %%  (1) some relativity is still contained
   %%  (2) relxillBB uses the proper alpha-disk Temperature Profile
   %%      and not the approximated one from diskbb
   if (diff > 0.01)
     status = EXIT_FAILURE;

%   ylog;xlog;
%   hplot(lo,hi,val);
%   ohplot(lo,hi,val_ref);

   
   return status;
}

