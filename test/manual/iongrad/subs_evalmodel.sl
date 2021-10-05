require("isisscripts");

load_xspec_local_models("build");



variable col_grad = ["#067bc2","#8b61c4","#d9258b","#e90c2a"];

define eval_mo(){
   variable lo,hi;
   (lo,hi) = log_grid(0.2,100,1000);
   variable val = eval_fun_keV(lo,hi)/(hi-lo)*lo*lo;
   return struct{lo=lo, hi=hi, val=val};
}

variable pl = xfig_plot_new();
variable pl_res = xfig_plot_new();
pl_res.xlabel("Energy [keV]");
pl_res.ylabel("Ratio");
pl.ylabel("$\nu F_\nu$"R);
variable colid=1;

define plot_mo(mo, val_ref){
   pl.world(min(mo.lo),max(mo.hi),10,350;loglog);
   pl_res.world(min(mo.lo),max(mo.hi),0.95,1.05;loglog);
   
   variable col = qualifier("color",colid);
   colid++;
   
   pl.plot(mo.lo,mo.val; color=col, loglog);
   pl_res.plot(mo.lo,mo.val/val_ref; color=col, loglog);

}


%%%%%%%%%% 
define relxill_read_iongrad_data(){
   variable rlo,rhi,lxi,logn;
   (rlo,rhi,lxi,logn) = readcol("__relxillOutput_iongrad.dat",1,2,3,4);
   return struct{rlo=rlo, rhi=rhi, lxi=lxi, logn=logn};
}


variable pl_lxi = xfig_plot_new();
pl_lxi.xlabel("radius [$R_\mathrm{g}$]"R);
pl_lxi.ylabel("$\log \xi$"R);

variable pl_ln = xfig_plot_new();
pl_ln.xlabel("radius [$R_\mathrm{g}$]"R);
pl_ln.ylabel("$\log N$"R);

define plot_grad(){

   () = eval_fun(1,2);
   variable str = relxill_read_iongrad_data();
   
   variable rmean = 0.5*(str.rlo+str.rhi);
   
          
   pl_lxi.world(min(str.rlo),max(str.rhi),0,5.0;xlog);
   pl_lxi.plot(rmean, str.lxi; color=qualifier("color", "black") );

   
   pl_ln.world(min(str.rlo),max(str.rhi),14,22;xlog);
   pl_ln.plot(rmean, str.logn; color=qualifier("color", "black") );

   
   variable rad_lxi_max = (11./9.)^2*str.rlo[0];
   pl_lxi.plot([rad_lxi_max, rad_lxi_max], [0,1]; world10, line=1);
   pl_ln.plot ([rad_lxi_max, rad_lxi_max], [0,1]; world10, line=1);

}