require("subs_evalmodel.sl");



define makeplot(h, fname){
      
   fit_fun("relxilllpCp");
   list_par;
   set_par("*iongrad_type",2);
   
   putenv("RELXILL_WRITE_OUTFILES=1");
   putenv("RELXILL_NUM_RZONES=50");   
   
   variable ii, n = length(h);
   variable leg = String_Type[n];
   variable col = col_grad;
   
   _for ii(0,n-1){
      set_par("*.h",h[ii],1,2,100);
      leg[ii] = sprintf("$h = %.1f$",h[ii]);
      plot_grad(;color=col[ii]);
            
   }
   
   pl_lxi.add_object(xfig_new_legend(leg,col,0,3,0.5),0.96,0.96,0.5,0.5;world0);
%C   pl.title(sprintf("$h=%.2f$",h));
   
   xfig_multiplot(pl_lxi,pl_ln;cols=1).render(fname);
   
}

variable h = [2.0,4.0,10,100];
variable fname = "plot_iongrad.pdf";
makeplot(h, fname);