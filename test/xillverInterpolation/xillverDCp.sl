require("isisscripts");

%use_localmodel("relxill");
load_xspec_local_models("./build/");


variable col_scheme = ["#03071e","#370617","#6a040f","#9d0208","#d00000","#dc2f02","#e85d04","#f48c06","#faa307","#ffba08"];
col_scheme = ["#007f5f","#2b9348","#55a630","#80b918","#aacc00","#bfd200","#d4d700","#dddf00","#eeef20","#ffff3f"];

variable lo,hi; 
(lo,hi) = log_grid(0.1,200,3000);
variable x = 0.5*(lo+hi);

define std_eval(){
   return eval_fun_keV(lo,hi)*x^2/(hi-lo);
}

%variable logn = [3.1,3.2,3.3,3.4,3.5,3.6,3.7]-3.0;

define std_model_plot(ff, pvals, pstring, plab){
   
   variable n = length(pvals);
   fit_fun(ff);

   variable ii;
   variable vals = Array_Type[n];
   
   variable pl = xfig_plot_new(13,10);
   pl.world(lo[0],hi[-1],qualifier("ymin",0.5),qualifier("ymax",3e3);loglog);
   
   variable lab = String_Type[n];
   
   _for ii(0,n-1){

      set_par(pstring,pvals[ii]);
      
      vals[ii] = std_eval();
      
      lab[ii] = sprintf(plab, pvals[ii]);
      
      pl.plot(x,vals[ii]; color=col_scheme[ii+1]);
   }
   
   pl.add_object( xfig_new_text(qualifier("title",ff)),0.5,0.98,0,0.5;world0);
   
   pl.xlabel("energy [keV]");
   pl.ylabel(`flux [$\nu F_\nu$] `);
   
   pl.add_object(xfig_new_legend(lab,col_scheme,0,3,0.5; labelsize="footnotesize"),
	      1.02,0.98,-0.5,0.5;world0);
      
   return pl;
}


variable logn = [15:16.2:#7];
variable pl1 = std_model_plot("xillverDCp", logn, "*.logN", `$\log N = %.2f$`);
variable logxi = [3:3.3:#7];
variable pl2 = std_model_plot("xillverDCp", logxi, "*.logxi", `$\log \xi = %.2f$`);

xfig_new_vbox_compound(pl1,pl2).render("plot_xillverDCp_detail.pdf");


#iffalse
variable logn = [15,15.5,16,16.5,17,17.5,18,18.5,19];
variable pl1 = std_model_plot("relxillDCp", logn, "*.logN", `$\log N = %.1f$`;ymin=10);
variable logxi = [0:4:#9];
variable pl2 = std_model_plot("relxillDCp", logxi, "*.logxi", `$\log \xi = %.2f$`;ymin=10);

xfig_new_vbox_compound(pl1,pl2).render("plot_relxillDCp.pdf");

