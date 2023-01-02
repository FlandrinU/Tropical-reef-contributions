
plan <- drake_plan(
  
  ##### set up output #####
  out_folder = dir.create("output"),
  out_folder_data = dir.create("output/data"),
  out_folder_plots = dir.create("output/plots"),

  ##### Bioenergetic model #####
  
  # load all parameters
  sptl = get_sptl(),
  kmax = get_kmax(),
  lw = get_lw(sptl),
  cnp = get_cnp(),
  cnpdiet = get_cnpdiet(),
  metpar = get_metpar(kmax, sptl),
  artr = get_artr(sptl),
  
  # combine parameters
  params = combine_params(sptl, kmax, lw, cnp, cnpdiet, metpar, artr),
  output_params = write.csv(params, "output/data/params_sst_glob.csv", row.names = FALSE),
  tfish_unique = get_tfish_unique(),
  
  # Run fishflux bioenergetic model
  cnpflux = run_fishflux(tfish_unique, params),
  output_cnpflux = write.csv(cnpflux, "output/data/flux_results.csv", row.names = FALSE),
  
  ##### Summarize transect level #####
  
  # load all transect data
  tfish = read.csv("data/tfish_corr.csv"),
  
  # summarize results per transect
  summary_transect = summarise_pertransect(tfish, cnpflux, params),
  
  contributions = get_contributions(tfish, cnpflux, params, summary_transect),
  contributions_site = get_contributions_site(tfish, cnpflux, params, summary_transect),
  sp_loc = unique(select(ungroup(contributions), bioregion, species)),
  
  # herbivory, piscivory 
  herb_pisc = get_herb_pisc(tfish, cnpflux, params, summary_transect),
  
  # combine summary data, sst 
  summary_transect_complete = combine_summary(summary_transect, herb_pisc),
  output_fluxsum = write.csv(summary_transect_complete, "output/data/flux_summary_transect.csv", row.names = FALSE),
  
  ##### models #####
  
  # site loc
  mod_mv_siteloc = fit_mvfun_siteloc(summary_transect_complete),

  # add biomass sst
  mod_mvfun_bm = fit_mvfun_bm(summary_transect_complete),

  # community analysis
  mod_mvfun_com = fit_mvfun_com(summary_transect_complete),
  mod_mvfun_com2 = fit_mvfun_com2(summary_transect_complete), #get absolute slopes

  # output tables 
  tab_mod_mv_siteloc = make_table_mod_mv_siteloc(mod_mv_siteloc),
  tab_mod_mvfun_bm = make_table_mod_mvfun_bm(mod_mvfun_bm),
  tab_mod_mvfun_com = make_table_mod_mvfun_com(mod_mvfun_com),

  ##### contribution analysis #####
  degree_dominance = get_dd_site(contributions_site),
  sp_importance = get_importance_site(contributions_site),
  freq_dominance = get_fd_site(sp_importance),
  
  mod_dd = fit_mod_dd(degree_dominance),
  mod_fd = fit_mod_fd(freq_dominance),

  # ##### FIGURES #####
  coords = get_coords_siteloc(summary_transect_complete),
  # # main
  fig1 = make_fig1(mod_mv_siteloc, mod_mvfun_bm, coords, summary_transect),
  fig2 = make_corplot(mod_mvfun_bm),
  
  fig3 = make_fig3(mod_mvfun_com),
  fig4 = make_fig4(degree_dominance, sp_importance, freq_dominance, mod_dd, mod_fd),
  # 
  # # Supplemental data figures
  SI_fig2 = make_pp_plots(mod_mv_siteloc, mod_mvfun_bm, mod_mvfun_com2)
)
