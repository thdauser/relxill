require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_model.sl");

load_relxill_dev_model();

_traceback=1;


variable params = Assoc_Type[Array_Type];

params["a"] = [0.0,0.2,0.4,0.5,0.6,0.7,0.9,0.95,0.99,0.998];
params["Rin"] = [-1,-1.1,-1.5,-2,-10];
params["Rout"] = [1000, 400, 200, 100];
params["Incl"] = [10:80:10];
params["logN"] = [15:19:1];
params["logxi"] = [1:3:#5];
params["z"] = [0:0.5:#5];
params["refl_frac"] = [0:5:1];



%%%% FUNCTIONS %%%%% 




variable pl = plotAllParameterCombinations(params);
variable combpl = xfig_multiplot(pl);
combpl.render(dir_plots+"plotModelParameterCombinations.pdf");
