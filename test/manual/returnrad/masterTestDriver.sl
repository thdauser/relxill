#!/usr/bin/env isis-script

_traceback=1;

%%% Tests %%% 
require("testRframeReturnradRefdata.sl");
require("testDiagnosePlots.sl");
require("testXillverPrimaryNormalization.sl");


if( testRframeReturnradRefdata() != EXIT_SUCCESS) exit(1);
if( testXillverPrimaryNormalization() != EXIT_SUCCESS) exit(1);
if( testDiagnosePlots() != EXIT_SUCCESS) exit(1);
