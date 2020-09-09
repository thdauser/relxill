variable ff_default = "relxillBBret";


define load_relxill_dev_model(){
   variable modlib = "model/librelxill.so";
   
   if (stat_file(modlib) == NULL){
      vmessage("\n **** error : local relxill %s/%s model not found  ; exiting ... **** \n ",
	     getcwd(), modlib);
      exit;
   }
   
   require("xspec");
   load_xspec_local_models(modlib);
}


define fit_fun_default(){ %{{{
   %% only works for single model components %%
   
   fit_fun(ff_default);
   
   variable p,pa = get_params();   
   variable df;
   foreach p(pa){
      df = eval(sprintf("%s_default(%i)",qualifier("ff0",ff_default),p.index-1));
      p.value=df.value;
      p.min=df.min;
      p.max=df.max;
      p.freeze=df.freeze;
   }
   set_params(pa);
}
%}}}


define get_ener_grid(){
   variable lo0, hi0;
   (lo0,hi0) = log_grid(0.1,50,3000);
   variable nuFnu = (0.5*(lo0+hi0)^2)/(hi0-lo0);
   return struct{lo=lo0,hi=hi0,val,nuFnu=nuFnu};
}


define get_std_plot(){
   
   variable ymin = 1e-3;
   variable ymax = 100;
   
   variable pl = xfig_plot_new();
   
   variable egrid = get_ener_grid();
   pl.world(egrid.lo[0],egrid.hi[-1],ymin,ymax;loglog);
   
   pl.xlabel("energy [keV]");
   pl.ylabel("flux [a.u.]");

   return pl;
}

define add_plot_info(pl,lab,key,col){
 
   pl.add_object(xfig_new_text(key), 0.5,0.98,0,0.5;world0);
   
   pl.add_object(xfig_new_legend(lab,col,0,3,0.5;fontsize="labelsize"),1.02,0.98,-0.5,0.5;world0);
      
}


define calc_single_evaluation(pl,egrid,key,value){
      
   set_par(sprintf("*.%s",key), value);
   
   egrid.val = eval_fun_keV(egrid.lo,egrid.hi); %*egrid.nuFnu;       

   return egrid;
   
}




define plot_single_param(pars,key){
   
   variable param_vals = pars[key];
   variable ii, n = length(param_vals);
   
   variable col = [CB_COLOR_SCHEME_NB,CB_COLOR_SCHEME_NB];
   
   variable pl = get_std_plot();
   
   variable lab = String_Type[n];
   variable keyTexSafe = strreplace(key,"_","-");
   
   fit_fun_default();
   variable egrid = get_ener_grid();
   _for ii(0,n-1){
           
      egrid = calc_single_evaluation(pl,egrid,key,param_vals[ii]);
      pl.plot(0.5*(egrid.lo+egrid.hi), egrid.val; color=col[ii]);

      lab[ii] = sprintf("%s = %.3e", keyTexSafe, param_vals[ii]);
   }
     
   add_plot_info(pl,lab,keyTexSafe,col);
   
   return pl;
}



define plotAllParameterCombinations(params){
   
   variable keys = assoc_get_keys(params);
   
   variable ii, n = length(keys);
   variable pl = Struct_Type[n];
   
   _for ii(0,n-1){      
      pl[ii] = plot_single_param(params,keys[ii]);
   }
   
   return pl;
}
