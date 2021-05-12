require("isisscripts");


variable dir = "../../../cmake-build-debug/";
variable lo,hi,val, val2;

(lo,hi,val) = readcol(dir+"test-debug-xillver-gshift-0.dat",1,2,3);

(lo,hi,val2) = readcol(dir+"test-debug-xillver-gshift-z.dat",1,2,3);

xlog;ylog;

hplot(lo,hi,val/(hi-lo));

ohplot(lo,hi,val2/(hi-lo));


#iffalse
hplot(lo,hi,val/(hi-lo));

ohplot(lo,hi,val2/(hi-lo));

ohplot(lo,hi,val/(hi-lo)*gshift);

ohplot(lo,hi,val/(hi-lo)*1.5);

ohplot(lo,hi,val/(hi-lo)/1.5);

hplot(lo,hi,val);

ohplot(lo,hi,val2);

fit_fun("powerlaw");

pl = eval_fun_keV(lo,hi);

list_par;

fit_fun("zpowerlw");

list_par;

pl = eval_fun_keV(lo,hi);

set_par(3,1/0.666-1);

list_par;

pl2 = eval_fun_keV(lo,hi);

ohplot(lo,hi,pl);

set_par(2);

pl2 = eval_fun_keV(lo,hi);

set_par(3,0);

pl = eval_fun_keV(lo,hi);

ohplot(lo,hi,pl);

set_par(2,2);

pl = eval_fun_keV(lo,hi);

set_par(3,0.5);

pl2 = eval_fun_keV(lo,hi);

ohplot(lo,hi,pl);

ohplot(lo,hi,pl*100);

ohplot(lo,hi,pl2*100);

1.66*1.5;
