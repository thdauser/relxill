%% energy shift for the returning radiation g = E_i / E_e = g_e / g_i
%% = ut_d(r_i) / ut_d(r_e)
define gfac_rr (ri, re, a, lambda){
   ut_d_lambda (ri, a, lambda) / ut_d_lambda (re, a, lambda);
}
