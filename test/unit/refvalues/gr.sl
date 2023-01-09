implements("gr");

define load_scripts_from_directory(scripts){   
   variable scr;   
   foreach scr(scripts){
      require(scr);
   }
}

variable gr_scripts = "gr/"+
  ["kerr.sl",
   "lamp_post.sl",
   "returning_radiation.sl",
   "disk.sl"   
  ];


load_scripts_from_directory(gr_scripts);



%!%+
%\function{GR Equations}
%\synopsis{contains many basic GR functions (Kerr, LP, ...)  }
%\usage{gr->ht_h(...)   (as loaded in namespace gr)
%\description
%   the following functions are loaded in the "gr" namespace
%
%   define ad_E(r,a)
%   define ad_Lz(r,a)
%   define ad_lambda(r,a)
%   define qrad_disk(r,rin)
%   define q2(del, h, a)
%   define v_eta (a,q2,eta)
%   define v_r (a,q2,r)
%   define v_phi(a,r)
%   define ut_d (r, a)
%   define ut_d_lambda (r, a, lambda)
%   define del_inc(q,r,a)
%   define cal_gamma(r,a)
%   define cal_a_gr(r,a)
%   define ut_h (h, a) 
%   define gi(h,r,a)
%   define ut_h_v(h, a, v)
%   define gravRedshift(height, a)
%   define gfac_rr (ri, re, a, lambda)   
%!%-
define gr(){
   help(_function_name()); 
}