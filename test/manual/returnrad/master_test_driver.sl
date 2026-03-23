#!/usr/bin/env isis-script

_traceback=1;

%%% Tests %%%
variable dir = "tests";
require("$dir/testRframeReturnradRefdata.sl"$);
require("$dir/testXillverPrimaryNormalization.sl"$);
require("$dir/testEmissivityBbody.sl"$);

variable status = EXIT_SUCCESS;
if( testRframeReturnradRefdata(;plot) != EXIT_SUCCESS) {
   vmessage(" *** error: test Restframe Returning Spectra failed \n\n");
   status = EXIT_FAILURE;
}
if( testXillverPrimaryNormalization() != EXIT_SUCCESS) {
   vmessage(" *** error: testing Xillver Primary Normalization failed \n\n");
   status = EXIT_FAILURE;
}
if( testEmissivityBbody() != EXIT_SUCCESS) {
   status = EXIT_FAILURE;
}


if (status==EXIT_FAILURE){
   vmessage("\n\n  ****** TESTS FAILED ******* \n\n");
   exit(1);
} 