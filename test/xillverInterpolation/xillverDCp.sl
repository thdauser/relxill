require("isisscripts");

require("subs.sl");


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

