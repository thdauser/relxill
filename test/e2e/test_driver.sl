#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

%%% OPTIONS %%%
variable verbose = 0;
#ifdef VERBOSE
verbose=1;
#endif


%%% "isis-script [-DSTABLE] test_driver.sl [tests/test_*.sl]"
%%% - setting the option -DSTABLE means that only the release models in
%%%  "lmodel_relxill_public.dat" are tested
%%% - if no script is given, all scripts with the above glob construct
%%%   are called 

variable TEST_DEVEL = 1;
#ifdef STABLE
TEST_DEVEL=0;
#endif


%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%
variable filename_tests = "tests/test_*.sl";

if (length(__argv)>1){
      filename_tests = __argv[1];
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


require("isisscripts");
require("subs_fitfunctions.sl");
require("load_test_setup.sl");

variable ALL_FF = get_implemented_fitfunctions(;dev=TEST_DEVEL);
Fit_Verbose=-1;

variable test_scripts = glob(filename_tests);
variable ntests = length(test_scripts);
variable retval = Int_Type[ntests];


define evaluate_single_test(fname_script){
   variable namespace = strchop(strchop(fname_script,'/',0)[-1], '.',0)[0];
   vmessage("[%s] %s ", namespace, fname_script);
   require(fname_script, namespace);

   variable retval = eval(sprintf("%s->runtest(ALL_FF);",namespace));
   
   if( (verbose || retval!=EXIT_SUCCESS)  ){
      eval(sprintf("if(msg_log!=NULL) {message(msg_log); }"), namespace);
   }
   
   return retval;
}



variable ii;
_for ii(0, ntests-1){

   retval[ii] = evaluate_single_test(test_scripts[ii]);
   
}




message( "\n==================== ");
vmessage(" Tests passed %i/%i ", nint(sum(retval)), ntests);
message( "==================== ");
variable passed = all(retval);
if (passed==0){
   message(" *** FAILURE *** ");
   exit(1);
} else {
   exit(0);
}
