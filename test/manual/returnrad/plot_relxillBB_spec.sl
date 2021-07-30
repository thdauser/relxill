require("isisscripts");


require("scripts/subs_returnrad.sl");
load_xspec_local_models("build/");

putenv("RELXILL_BBRET_NOREFL=1");

fit_fun("relxillBB");
variable spin = 0.0;
variable Tin = 1.0;

set_par("*.kTbb",Tin);

list_par;
set_par("*.a",spin);

variable lo, hi;
(lo, hi) = log_grid(0.5,20,1000);

set_par("*.boost",-1,0,-10,10);
variable refl = eval_fun_keV(lo,hi)/ (hi-lo)*lo;

set_par("*.boost",0);
variable prim = eval_fun_keV(lo,hi)/ (hi-lo)*lo;

set_par("*.boost",1);
variable comb = eval_fun_keV(lo,hi)/ (hi-lo)*lo;

xlog;ylog;
hplot(lo,hi,prim);
ohplot(lo,hi,refl);
ohplot(lo,hi,comb);



fit_fun("diskbb");
set_par("*.Tin",1.0);
variable spec_disk = eval_fun_keV(lo,hi)/ (hi-lo)*lo;
ohplot(lo,hi,spec_disk*sum(prim)/sum(spec_disk));

%variable spec_kerrbb = eval_kerrbb_model(lo,hi,spin);
%ohplot(lo,hi,spec_kerrbb*sum(prim)/sum(spec_kerrbb));
