require("subs_evalmodel.sl");



define makeplot(h, fname){
      
   fit_fun("relxilllpCp");
   list_par;
   set_par("*iongrad_type",2);
   set_par("*.h",h);
   set_par("*.logN",19);
   
   putenv("RELXILL_WRITE_OUTFILES=1");
   
   putenv("RELXILL_NUM_RZONES=50");   
   variable mo_ref = eval_mo();   
   
   
   variable ii, nzones = [50,30,20,15,10];
   variable n = length(nzones);

   variable leg = String_Type[n];
   variable col = ["black",col_grad];
   
   _for ii(0,n-1){
      putenv(sprintf("RELXILL_NUM_RZONES=%i",nzones[ii]));
      leg[ii] = sprintf("N = %i",nzones[ii]);
      plot_mo(eval_mo(), mo_ref.val;color=col[ii]);
      
      
   }
   
   pl.add_object(xfig_new_legend(leg,col,0,3,0.5),0.96,0.96,0.5,0.5;world0);
   pl.title(sprintf("$h=%.2f$",h));
   
   xfig_multiplot(pl,pl_res;cols=1).render(fname);
   
}

variable h = 3.0;
variable fname = "plot_iongrad_numzones.pdf";
makeplot(h, fname);