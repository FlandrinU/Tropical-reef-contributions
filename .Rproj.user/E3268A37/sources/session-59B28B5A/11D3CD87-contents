
make_fig1 <- function(mod_mv_siteloc, mod_mvfun_bm, coords, summary_transect){
  
  nsiteloc <- summary_transect %>%
    group_by(locality) %>%
    mutate(nsite = length(unique(sites))) %>%
    select(locality, sites, nsite) %>%
    unique() 
  
  ls1 <- nsiteloc[nsiteloc$nsite==1, ]$locality
  ls2 <- nsiteloc[nsiteloc$nsite==1, ]$sites
  
  re1 <- ranef(mod_mv_siteloc)
  re_loc1 <- re1$locality[,1,] %>%
    as.data.frame()
  re_site1 <- re1$sites[,1,] %>%
    as.data.frame()
  re_loc1[ls1,] <- re_loc1[ls1,] + re_site1[ls2,]
  
  colnames(re_loc1) <- c("r1_loc_Fn", "r1_loc_Fp", "r1_loc_Gc", "r1_loc_Iherb", "r1_loc_Ipisc")
  re_loc1 <- re_loc1 %>%
    rownames_to_column("locality") 
  
  re_loc1 <- re_loc1 %>%
    group_by() %>%
    mutate(r1 = normalize(r1_loc_Fn),
           r2 = normalize(r1_loc_Fp),
           r3 = normalize(r1_loc_Gc),
           r4 = normalize(r1_loc_Iherb),
           r5 = normalize(r1_loc_Ipisc)) %>%
    ungroup()
  
  re_loc1$r_multi <- sapply(1:98, function(i){
    re_loc1[i,] %>%
      select(r1, r2, r3, r4, r5) %>%
      simplify() %>%
      geomean()
  })


  coords <- coords %>%
    select(- sites, -lat, -lon, -mean) %>%
    unique()
  
  data <- re_loc1 %>%
    left_join(coords) %>%
    unique()
  
  library("rnaturalearth")
  library("rnaturalearthdata")
  library(rgeos)
  library(sf)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  class(world)
  
  
  
  theme_worldmap <- function(){
    theme_bw() +
      theme(legend.position = "none", 
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA, color = "grey", size = 1),
            axis.title = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6),
            plot.title = element_text(size = 6, hjust = 0, face = "bold"), 
            plot.margin = unit(c(0.002,0.002,0.002,0.002), units = , "cm"),
            panel.spacing = unit(0, "cm")
      ) 
  }
  
  col <- fish(n = 5, option = "Callanthias_australis")
  
  get_scale <- function(f){
    quant <- quantile(f, c(0,0.25,0.75,1), na.rm = TRUE)
    quant[1] <- quant[1] - 0.5   # adjust lower bound
    factor(cut(f, breaks = quant), labels = c( "low", "medium", "high"))
  }
  
  get_scores <- function(proc){
    probs <- seq(0,1,0.01)
    quantiles <- quantile(proc, prob = probs, na.rm = TRUE)
    quant <- as.numeric(factor(findInterval(proc, quantiles)))/100
    return(quant)
  }
  

  getcol <- colorRampPalette(c("white", col[1]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","white"))
  pal1 <- getcol(n = 3)
  pal <- c(pal1[1], pal2[2:3])
  pal <- c("black", "grey60", col[1])
  
  
  
  a <-
    ggplot(data) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon_loc, y = lat_loc, color = get_scale(r1_loc_Fn), 
                   size = get_scores(r1_loc_Fn)),  alpha = 0.9,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon_loc, y = lat_loc),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(data, r1_loc_Fn > quantile(data$r1_loc_Fn, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, 
                  label = as.character(expression(paste("N excretion (gN ", m^{-2}, day^{-1}, ")")))), 
              size = 2, hjust = 0, parse = T) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm"))
  a
  
  
  getcol <- colorRampPalette(c("white", col[2]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","white"))
  pal1 <- getcol(n = 3)
  pal <- c("black", "grey60", col[2])
  
  
  b <-
    ggplot(data) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon_loc, y = lat_loc, color = get_scale(r1_loc_Fp), 
                   size = get_scores(r1_loc_Fp)),  alpha = 0.9,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon_loc, y = lat_loc),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(data, r1_loc_Fp > quantile(data$r1_loc_Fp, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, 
                  label = as.character(expression(paste("P excretion (gP ", m^{-2}, day^{-1}, ")")))), 
              size = 2, hjust = 0, parse = T) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.0001,0.0001,0.0001,0.0001), units = , "cm"))
  b
  
  
  
  getcol <- colorRampPalette(c("white", col[3]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","white"))
  pal1 <- getcol(n = 3)
  pal <- c("black", "grey60", col[3])
  
  c <-
    ggplot(data) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon_loc, y = lat_loc, color = get_scale(r1_loc_Gc), 
                   size = get_scores(r1_loc_Gc)),  alpha = 0.9,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon_loc, y = lat_loc),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(data, r1_loc_Gc > quantile(data$r1_loc_Gc, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, 
                  label = as.character(expression(paste("Production (gC ", m^{-2}, day^{-1}, ")")))), 
              size = 2, hjust = 0, parse = T) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.9))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.00,0.00,0.00,0.00), units = , "cm"))
  c
  
  
  
  
  
  
  getcol <- colorRampPalette(c("lightgrey", col[4]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","lightgrey"))
  pal1 <- getcol(n = 3)
  pal <- c(pal1[1:2], pal2[2:3])
  pal <- c("black", "grey60", col[4])
  
  
  d <-
    ggplot(data) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon_loc, y = lat_loc, color = get_scale(r1_loc_Iherb), 
                   size = get_scores(r1_loc_Iherb)),  alpha = 0.9,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon_loc, y = lat_loc),  alpha = 0.9,
               shape = 1, size = 4,
               data = filter(data, r1_loc_Iherb > quantile(data$r1_loc_Iherb, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, label = 
                    as.character(expression(paste("Herbivory (gC ", m^{-2}, day^{-1}, ")")))), 
              size = 2, hjust = 0, parse = T) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.00,0.00,0.00,0.00), units = , "cm"))
  d
  
  
  
  getcol <- colorRampPalette(c("lightgrey", col[5]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","lightgrey"))
  pal1 <- getcol(n = 3)
  pal <- c(pal1[1:2], pal2[2:3])
  
  pal <- c("black", "grey60", col[5])
  
  e <-
    ggplot(data) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon_loc, y = lat_loc, color = get_scale(r1_loc_Ipisc), 
                   size = get_scores(r1_loc_Ipisc)),  alpha = 0.9,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon_loc, y = lat_loc),  alpha = 0.9,
               shape = 1, size = 4,
               data = filter(data, r1_loc_Ipisc > quantile(data$r1_loc_Ipisc, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "",
                       drop = TRUE, na.translate = F) +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, 
                  label = as.character(expression(paste("Piscivory (gC ", m^{-2}, day^{-1}, ")")))), 
              size = 2, hjust = 0, parse = T) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_worldmap()+ 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.0,0.00,0.00,0.00), units = , "cm"))
  e
  

  
  pal <- c("black", "grey60", "deepskyblue1")
  
  f <-
    ggplot(data) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon_loc, y = lat_loc, color = get_scale(r_multi), 
                   size = get_scores(r_multi)),  alpha = 0.9,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon_loc, y = lat_loc),  alpha = 0.9,
               shape = 1, size = 4,
               data = filter(data, r_multi > quantile(data$r_multi, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "",
                       drop = TRUE, na.translate = F) +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, label = "Multifunctionality"), size = 2, hjust = 0) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_worldmap()+ 
    theme(legend.position = "none", 
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.0,0.00,0.00,0.00), units = , "cm"))
  f
  
  #### relationship with biomass 
  
  
  nd <- data.frame(biomass_tot = seq(50, 500, 5),
                   mean = 27,
                   locality = NA,
                   sites = NA)
  
  pred <- predict(mod_mvfun_bm, newdata = nd)
  
  ef <- predict(mod_mvfun_bm, newdata = nd, summary = F, nsamples = 1000)

  ef[,,1] <- normalize(exp(ef[,,1]))
  ef[,,2] <- normalize(exp(ef[,,2])) 
  ef[,,3] <- normalize(exp(ef[,,3])) 
  ef[,,4] <- normalize(exp(ef[,,4])) 
  ef[,,5] <- normalize(exp(ef[,,5])) 
  
  mf <- apply(ef, 1:2, geomean)

  mf_m <- apply(mf, 2, mean) 
  mf_lb <- apply(mf, 2, quantile, 0.025) 
  mf_ub <- apply(mf, 2, quantile, 0.975) 
  

  
  a1 <- 
    ggplot() +
    geom_ribbon(aes(x = nd$biomass_tot, ymin = exp(pred[,3,1]), ymax = exp(pred[,4,1])), 
                alpha = 0.5, fill = col[1]) +
    geom_smooth(aes(x = nd$biomass_tot, y = exp(pred[,1,1])), alpha = 0.5, color = col[1]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Biomass (g/m2)", y = "N excretion (gN/m2/day)") +
    theme(axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  b1 <-
    ggplot() +
    geom_ribbon(aes(x = nd$biomass_tot, ymin = exp(pred[,3,2]), ymax = exp(pred[,4,2])), 
                alpha = 0.5, fill = col[2]) +
    geom_smooth(aes(x = nd$biomass_tot, y = exp(pred[,1,2])), alpha = 0.5, color = col[2]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Biomass (g/m2)", y = "P excretion (gP/m2/day)") +
    theme(axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  c1 <- 
    ggplot() +
    geom_ribbon(aes(x = nd$biomass_tot, ymin = exp(pred[,3,3]), ymax = exp(pred[,4,3])), 
                alpha = 0.5, fill = col[3]) +
    geom_smooth(aes(x = nd$biomass_tot, y = exp(pred[,1,3])), alpha = 0.5, color = col[3]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Biomass (g/m2)", y = "Production (gC/m2/day)") +
    theme(axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  d1 <- 
    ggplot() +
    geom_ribbon(aes(x = nd$biomass_tot, ymin = exp(pred[,3,4]), ymax = exp(pred[,4,4])), 
                alpha = 0.5, fill = col[4]) +
    geom_smooth(aes(x = nd$biomass_tot, y = exp(pred[,1,4])), alpha = 0.5, color = col[4]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Biomass (g/m2)", y = "Herbivory (gC/m2/day)") +
    theme(axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  e1 <- 
    ggplot() +
    geom_ribbon(aes(x = nd$biomass_tot, ymin = exp(pred[,3,5]), ymax = exp(pred[,4,5])), 
                alpha = 0.5, fill = col[5]) +
    geom_smooth(aes(x = nd$biomass_tot, y = exp(pred[,1,5])), alpha = 0.5, color = col[5]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Biomass (g/m2)", y = "Piscivory (gC/m2/day)") +
    theme(axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  
  f1 <- 
    ggplot() +
    geom_ribbon(aes(x = nd$biomass_tot, ymin = mf_lb, ymax = mf_ub), 
                alpha = 0.5, fill = "deepskyblue1") +
    geom_smooth(aes(x = nd$biomass_tot, y = mf_m), alpha = 0.5, color = "deepskyblue1") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = (expression(paste("Biomass (g ", m^{-2}, ")"))), 
         y = "Multifunction", parse = T) +
    theme(axis.title.y = element_blank() )
  
  
  
  multimap <-
    a + a1 + b + b1 + c + c1 + d + d1 + e + e1 + f + f1 +
    plot_layout(ncol = 2, widths = c(4,1))# + plot_annotation(tag_levels = "a")
  
  #multimap
  ggsave("output/plots/fig1_multimap.png", multimap, width = 180, height = 210, units = "mm")
  ggsave("output/plots/fig1_multimap.pdf", multimap, width = 180, height = 210, units = "mm")
}



make_fig3 <- function(mod_mvfun_com){
  
  fe <- brms::fixef(mod_mvfun_com) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("dep", "var", "extra"), sep = "_") %>%
    dplyr::mutate(var = paste0(var, extra)) %>%
    dplyr::filter(!var == "InterceptNA") %>%
    filter(!var == "standardlogbiomasstot") %>%
    filter(!var == "standardmeanNA") %>%
    mutate(var = fct_relevel(var,
                             "standardimmq1",
                             "standardimmm",
                             "standardimmq3",
                             "standardsizeq1",
                             "standardsizem",
                             "standardsizeq3",
                             "standardtrophq1",
                             "standardtrophm",
                             "standardtrophq3",
                             "standardnspecNA"
    ))
  
  unique(fe$var)
  plot <-
    ggplot(fe) +
    geom_hline(yintercept = 0) +
    geom_linerange(aes(x = var, ymin = Q2.5, ymax = Q97.5, color = dep),
                   position = position_dodge(.7), size = 1, linetype = 1) +
    geom_point(aes(x = var, y = Estimate, color = dep), position = position_dodge(.7), size = 3) +
    coord_flip() +
    theme_bw() +
    labs(x = "", color = "") +
    theme(legend.position = "bottom", axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12),
          axis.line.x = element_line(), axis.ticks.x = element_line(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()
    )  
  
  slopes <-
    ggplot(fe) +
    geom_vline(xintercept = seq(1, 10 ,1) + .5, color = "lightgrey",
               size = 0.5, linetype = 1, alpha = 0.7) +
    geom_hline(yintercept = 0, size = 1, color = "black") +
    geom_linerange(aes(x = var, ymin = Q2.5, ymax = Q97.5, color = dep),
                   position = position_dodge(.7), size = 1, linetype = 1) +
    geom_point(aes(x = var, y = Estimate, color = dep), position = position_dodge(.7), size = 2) +
    coord_flip() +
   # scale_y_continuous(lim = c(-0.22, 0.22)) +
    theme_minimal() +
    scale_color_fish_d(option = "Callanthias_australis",
                       labels = c("N excretion", "P excretion", "Production", "Herbivory", "Piscivory"),
                       name = "")  +
    labs(y = "Effect", x = "") +
    scale_x_discrete(
      labels = rev(c( "Richness", "Trophic level (upper)", "Trophic level (median)", "Trophic level (lower)",
                      "Size (upper)", "Size (median)", "Size (lower)",
                      "Immaturity (upper)", "Immaturity (median)", "Immaturity (lower)"
      ))
    ) +
    
    theme(legend.position = "bottom", axis.text = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 12), legend.text = element_text(size = 10),
          axis.line.x = element_line(), axis.ticks.x = element_line(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank() 
    ) + guides(color=guide_legend(nrow=2, byrow=TRUE))

  slopes
  ggsave("output/plots/fig3_com_slopes.png", slopes, width = 180, height = 200, units ="mm") 
  ggsave("output/plots/fig3_com_slopes.pdf", slopes, width = 180, height = 200, units ="mm")
  
  
}


make_fig4 <- function(degree_dominance, sp_importance, freq_dominance, mod_dd, mod_fd){

  # simulating communities
  rank_area <- function(y) {
    n <- length(y)
    a <- y[2:n]
    b <- y[1:(n-1)]
    area <- sum(a + b)/2
    return(area)
  }
  ks_area <- function(data){
    
    if (length(data) == 0){
      return(NA)
    }else if (length(data) == 1){
      return(1)
    }else{
      data <- data.frame(data = data)
      drank <- data %>% 
        mutate(rank = dense_rank(desc(data))) %>%
        dplyr::arrange(rank) %>%
        mutate(cumsum = cumsum(data))
      
      A <- rank_area(drank$cumsum)
      
      r <- nrow(data) # species richness
      A_max <- 1 * (r - 1) # maximum if first species performs 100% (rectangle)
      A_min <- (r^2-1)/(2*r) # minimum if all sp contribute equally (trapezium)
      
      dd <- (A - A_min)/(A_max - A_min)
      
      return(dd)
    }
  }
  
  list <- parallel::mclapply(3:277, function(j){
    te <- (hitandrun::simplex.sample(j, 500, sort=FALSE))[[1]]
    dd <- (apply(te, 1, ks_area))
    data.frame(nspec = j, dd = dd)
  }, mc.cores = 50)
  
  vec <- plyr::ldply(list) %>%
    group_by(nspec) %>%
    dplyr::summarize(dd = mean(dd))
  
  summary(vec$dd)
  
  plot(vec$nspec, vec$dd)
  
  hist(vec$dd)
  
  nd <- data.frame(data=NA)
  
  pr1 <- fitted(mod_dd[[1]], newdata = nd, nsamples = 1000)
  pr2 <- fitted(mod_dd[[2]], newdata = nd, nsamples = 1000)
  pr3 <- fitted(mod_dd[[3]], newdata = nd, nsamples = 1000)
  pr4 <- fitted(mod_dd[[4]], newdata = nd, nsamples = 1000)
  pr5 <- fitted(mod_dd[[5]], newdata = nd, nsamples = 1000)
  
  ###### degree_dominance ######
 
  inverse_logit <- function(x){
    exp(x)/(1+exp(x))
  }
  

  cols <- fish(option = "Callanthias_australis", n = 5)
  
  plot1 <-
    ggplot(degree_dominance) + 
    geom_jitter(aes(y = "e", x = dd_Fn), color = cols[1],
                alpha = 0.5, size = 0.2) +
    geom_jitter(aes(y = "d", x = dd_Fp), color = cols[2],
                alpha = 0.5, size = 0.2) +
    geom_jitter(aes(y = "c", x = dd_Gc), color = cols[3],
               alpha = 0.5, size = 0.2) +
    geom_jitter(aes(y = "b", x = dd_I_herb), color = cols[4],
                 alpha = 0.5, size = 0.2, width = 0.01) +
    geom_jitter(aes(y = "a", x = dd_I_pisc), color = cols[5],
                 alpha = 0.5, size = 0.2, width = 0.01) +
    geom_pointrange(aes(y = "e", x = pr1[1], xmin = pr1[3], xmax = pr1[4])) +
    geom_pointrange(aes(y = "d", x = pr2[1], xmin = pr2[3], xmax = pr2[4])) + 
    geom_pointrange(aes(y = "c", x = pr3[1], xmin = pr3[3], xmax = pr3[4])) +
    geom_pointrange(aes(y = "b", x = pr4[1], xmin = pr4[3], xmax = pr4[4])) +
    geom_pointrange(aes(y = "a", x = pr5[1], xmin = pr5[3], xmax = pr5[4])) +
    geom_vline(xintercept = mean(vec$dd), linetype = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_discrete(labels = c("Piscivory", 
                                "Herbivory", 
                                "Production", 
                                "P excretion", 
                                "N excretion")) +
    labs(x = "Degree of dominance (communities)", y = "") +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line())
  
  spi_sp_long <- freq_dominance %>%
    pivot_longer(c(Fn_i, Fp_i, Gc_i, I_herb_i, I_pisc_i)) %>%
    drop_na() %>%
    filter(occ>4, value > 0) 
  
  a <- spi_sp_long %>% group_by(name) %>%
    summarize(n = n())
  
  b <- freq_dominance %>%
    pivot_longer(c(Fn_i, Fp_i, Gc_i, I_herb_i, I_pisc_i)) %>%
    drop_na() %>%
    filter(occ>4) %>% 
    group_by(name) %>%
    summarize(ntot = n())
  
  ab <- left_join(a,b) %>%
    mutate(prop = n/ntot)
  
  # get value of prop species important sometimes for a function
  freq_dominance %>%
    mutate(imp = Fn_i>0|Fp_i>0|Gc_i>0|I_pisc_i>0|I_herb_i>0) %>%
    group_by() %>%
    summarize(prop = sum(imp, na.rm = T)/n())
  # 49%
  
  
  plot2 <- 
    ggplot(ab) +
    geom_col(aes(x = prop, y = reorder(name, desc(name)), fill = name),
             alpha = 0.9) +
    geom_text(aes(x = prop - 0.1, y = reorder(name, desc(name)), 
                  label = paste0(round(100*prop), " %")),
              size = 4) +
    scale_fill_fish_d(option = "Callanthias_australis") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_discrete(labels = c("Piscivory", 
                                "Herbivory", 
                                "Production", 
                                "P excretion", 
                                "N excretion")) +
    scale_color_fish_d(option = "Callanthias_australis") +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "Proportion of species being dominant") +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line())
  
  
  pre1 <- fitted(mod_fd[[1]], newdata = data.frame(name = NA))
  pre2 <- fitted(mod_fd[[2]], newdata = data.frame(name = NA))
  pre3 <- fitted(mod_fd[[3]], newdata = data.frame(name = NA))
  pre4 <- fitted(mod_fd[[4]], newdata = data.frame(name = NA))
  pre5 <- fitted(mod_fd[[5]], newdata = data.frame(name = NA))
  
  plot3 <- 
    ggplot() + 
    geom_jitter(aes(y = reorder(name, desc(name)), x = value, color = name),
                data = spi_sp_long, alpha = 0.5, size = 0.2) +
    geom_pointrange(aes(y = 5, x = pre1[1,1], xmin = pre1[1,3], xmax = pre1[1,4])) +
    geom_pointrange(aes(y = 4, x = pre2[1,1], xmin = pre2[1,3], xmax = pre2[1,4])) +
    geom_pointrange(aes(y = 3, x = pre3[1,1], xmin = pre3[1,3], xmax = pre3[1,4])) +
    geom_pointrange(aes(y = 2, x = pre4[1,1], xmin = pre4[1,3], xmax = pre4[1,4])) +
    geom_pointrange(aes(y = 1, x = pre5[1,1], xmin = pre5[1,3], xmax = pre5[1,4])) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_discrete(labels = c("Piscivory", 
                                "Herbivory", 
                                "Production", 
                                "P excretion", 
                                "N excretion")) +
    scale_color_fish_d(option = "Callanthias_australis") +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "Frequency of dominance per species") +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line())
  
  plot <- plot1 + plot2 + plot3 + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'), axis.title.x = element_text(size = 9))
  
  ggsave("output/plots/fig4_dd.png", plot, width = 180, height = 90, units = "mm")
  ggsave("output/plots/fig4_dd.pdf", plot, width = 180, height = 90, units = "mm")
  
  
}


########### Figures supplemental ###########

make_pp_plots <- function(mod_mv_siteloc, mod_mvfun_bm, mod_mvfun_com2){
  
  a1 <- pp_check(mod_mv_siteloc, resp = "logFn", nsamples = 50) +
    labs(x = expression(paste("N excretion (gN ", m^{-1}, day^{-1}, ")")))
  a2 <- pp_check(mod_mv_siteloc, resp = "logFp", nsamples = 50) +
    labs(x = expression(paste("P excretion (gP ", m^{-1}, day^{-1}, ")")))
  a3 <- pp_check(mod_mv_siteloc, resp = "logGc", nsamples = 50) +
    labs(x = expression(paste("Production (gC ", m^{-1}, day^{-1}, ")")))
  a4 <- pp_check(mod_mv_siteloc, resp = "logIherb", nsamples = 50) +
    labs(x = expression(paste("Herbivory (gC ", m^{-1}, day^{-1}, ")")))
  a5 <- pp_check(mod_mv_siteloc, resp = "logIpisc", nsamples = 50) +
    labs(x = expression(paste("Piscivory (gC ", m^{-1}, day^{-1}, ")")))
  
  b1 <- pp_check(mod_mvfun_bm, resp = "logFn", nsamples = 50) +
    labs(x = expression(paste("N excretion (gN ", m^{-1}, day^{-1}, ")")))
  b2 <- pp_check(mod_mvfun_bm, resp = "logFp", nsamples = 50) +
    labs(x = expression(paste("P excretion (gP ", m^{-1}, day^{-1}, ")")))
  b3 <- pp_check(mod_mvfun_bm, resp = "logGc", nsamples = 50) +
    labs(x = expression(paste("Production (gC ", m^{-1}, day^{-1}, ")")))
  b4 <- pp_check(mod_mvfun_bm, resp = "logIherb", nsamples = 50) +
    labs(x = expression(paste("Herbivory (gC ", m^{-1}, day^{-1}, ")")))
  b5 <- pp_check(mod_mvfun_bm, resp = "logIpisc", nsamples = 50) +
    labs(x = expression(paste("Piscivory (gC ", m^{-1}, day^{-1}, ")")))
  

  c1 <- pp_check(mod_mvfun_com2, resp = "logFn", nsamples = 50) +
    labs(x = expression(paste("N excretion (gN ", m^{-1}, day^{-1}, ")")))
  c2 <- pp_check(mod_mvfun_com2, resp = "logFp", nsamples = 50) +
    labs(x = expression(paste("P excretion (gP ", m^{-1}, day^{-1}, ")")))
  c3 <- pp_check(mod_mvfun_com2, resp = "logGc", nsamples = 50) +
    labs(x = expression(paste("Production (gC ", m^{-1}, day^{-1}, ")")))
  c4 <- pp_check(mod_mvfun_com2, resp = "logIherb", nsamples = 50) +
    labs(x = expression(paste("Herbivory (gC ", m^{-1}, day^{-1}, ")")))
  c5 <- pp_check(mod_mvfun_com2, resp = "logIpisc", nsamples = 50) +
    labs(x = expression(paste("Piscivory (gC ", m^{-1}, day^{-1}, ")")))
  
  pp_plots <- a1 + a2 + a3 + a4 + a5 + b1 + b2 + b3 + b4 + b5 + c1 + c2 + c3 + c4 + c5 + 
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'a', tag_suffix = ")") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave("output/plots/pp_plots.png", pp_plots, width = 18, height = 12)
  ggsave("output/plots/pp_plots.pdf", pp_plots, width = 18, height = 12)
  
}

####### extra analysis diets #######
alt_diet <- function(tfish, cnpflux){
  
  diet <- read.csv("data/extrapolation_trophic_guilds.csv") 
  
  diet_experts <- read.csv("data/original_experts_classification.csv")%>%
    select(species = Genus_and_species_, Gaspar, Siquiera)
  
  diets <- diet %>%
    select(family, species, p2_m, p4_m, trophic_guild_predicted) %>% 
    left_join(diet_experts) %>%
    mutate(herb = trophic_guild_predicted == 2,
           herb_gas = Gaspar == "HD",
           herb_siq = Siquiera == "HD",
           herb_t1 = trophic_guild_predicted == 2 & p2_m > 0.5,
           pisc = trophic_guild_predicted == 4,
           pisc_gas = Gaspar == "FC",
           pisc_siq = Siquiera == "GC",
           pisc_t1 = trophic_guild_predicted == 4 & p4_m > 0.5
    )
  
  
  tfish <- left_join(tfish, unique(dplyr::select(cnpflux, species, v_m, size_cm, (ends_with("_median"))))) %>% 
    left_join(diets) %>%
    mutate(Ic = Ic_median * abun) %>%
    mutate(I = Ic/area) 
  
  
  tfs <- dplyr::group_by(tfish, region, locality, sites, area, transect_id) %>% 
    dplyr::summarise(
      I_herb = sum(I[herb == TRUE], na.rm = TRUE),
      I_pisc = sum(I[pisc == TRUE], na.rm = TRUE),
      I_herb_t1 = sum(I[herb_t1 == TRUE], na.rm = TRUE),
      I_pisc_t1 = sum(I[pisc_t1 == TRUE], na.rm = TRUE),
      I_herb_gas= sum(I[herb_gas == TRUE], na.rm = TRUE),
      I_pisc_gas = sum(I[pisc_gas == TRUE], na.rm = TRUE),
      I_herb_siq = sum(I[herb_siq == TRUE], na.rm = TRUE),
      I_pisc_siq = sum(I[pisc_siq == TRUE], na.rm = TRUE)
    ) %>% ungroup() %>%
    mutate(d_herb1 =  (I_herb_t1 - I_herb),
           d_pisc1 =  (I_pisc_t1 - I_pisc),
           d_herb_gas = (I_herb_gas - I_herb),
           d_pisc_gas = (I_pisc_gas - I_pisc),
           d_herb_siq = (I_herb_siq - I_herb),
           d_pisc_siq = (I_pisc_siq - I_pisc))
  
  b1 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_herb_gas), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-5, - 2.5, 0, 2.5, 5), limits = c(-5, 5)) +
    labs(y = "Percent of total", x = "Difference herbivory rate", title = "C) Comparison Mouillot") +
    theme_bw()
  b2 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_pisc_gas), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-1, - 0.5, 0, 0.5, 1), limits = c(-1, 1)) +
    labs(y = "Percent of total", x = "Difference piscivory rate", title = "D) Comparison Mouillot") +
    theme_bw()
  c1 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_herb_siq), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-5, - 2.5, 0, 2.5, 5), limits = c(-5, 5)) +
    labs(y = "Percent of total", x = "Difference herbivory rate", title = "E) Comparison Siqueira") +
    theme_bw()
  c2 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_pisc_siq), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-1, - 0.5, 0, 0.5, 1), limits = c(-1, 1)) +
    labs(y = "Percent of total", x = "Difference piscivory rate", title = "F) Comparison Siqueira") +
    theme_bw()
  
  p <-  b1 + b2 + c1 + c2 +plot_layout(ncol = 2)
  
  ggsave("output/plots/Supplemental_herbivory_piscivory_alternative_classification.pdf", 
         plot = p, 
         height = 6, width = 8)
  
}



##### correlations #####

make_corplot <- function(mod_mvfun_bm){
  
  sum <- summary(mod_mvfun_bm)
  rescor <- sum$rescor_pars
  loc <- sum$random$locality[-c(1:5),]
  site <- sum$random$sites[-c(1:5),]
  
  col <- fish(n = 10, option = "Ostracion_whitleyi", begin = 0, end = 1)
  
  # nd <- data.frame(biomass_tot = 100, mean = 27, locality = NA, sites = NA)
  # 
  # pred <- fitted(mod_mvfun_bm, newdata = nd, summary = FALSE, nsamples = 1000)

  pred <- residuals(mod_mvfun_bm)
  
  
  p1 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,2]), x = (pred[,1,1])),
               color = col[1], alpha = 0.2, size = 0.5) +
    theme_bw() +
    labs(x = "N excretion (gN/m2/day)", y = "P excretion (gP/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))
  p2 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,3]), x = (pred[,1,1])),
               color = col[2], alpha = 0.2, size = 0.5) +
    theme_bw() +
    labs(x = "N excretion (gN/m2/day)", y = "Production (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))
  p3 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,4]), x = (pred[,1,1])),
               color = col[3], alpha = 0.2, size = 0.5) +
    theme_bw() +
    labs(x = "N excretion (gN/m2/day)", y = "Herbivory (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))
  p4 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,5]), x = (pred[,1,1])),
               color = col[4], alpha = 0.2, size = 0.5) +
    theme_bw() +
    labs(x = "N excretion (gN/m2/day)", y = "Piscivory (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))
  p5 <-
  ggplot() +
    geom_point(aes(y = (pred[,1,3]), x = (pred[,1,2])),
               color = col[5], alpha = 0.2, size = 0.5) +
    theme_bw() +
    labs(x = "P excretion (gP/m2/day)", y = "Production (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(size = 8))
  p6 <-
  ggplot() +
    geom_point(aes(y = (pred[,1,4]), x = (pred[,1,2])),
               color = col[6], alpha = 0.2, size = 0.5)  +
    theme_bw() +
    labs(x = "P excretion (gP/m2/day)", y = "Herbivory (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))
  p7 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,5]), x = (pred[,1,2])),
               color = col[6], alpha = 0.2, size = 0.5)  +
    theme_bw() +
    labs(x = "P excretion (gP/m2/day)", y = "Piscivory (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(), 
          axis.title = element_text(size = 8))
  p8 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,4]), x = (pred[,1,3])),
               color = col[8], alpha = 0.2, size = 0.5)  +
    theme_bw() +
    labs(x = "Production (gC/m2/day)", y = "Herbivory (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))
  p9 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,5]), x = (pred[,1,3])),
               color = col[9], alpha = 0.2, size = 0.5)  +
    theme_bw() +
    labs(x = "Production (gC/m2/day)", y = "Piscivory (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(), 
          axis.title = element_text(size = 8))
  p10 <-
    ggplot() +
    geom_point(aes(y = (pred[,1,5]), x = (pred[,1,4])),
               color = col[10], alpha = 0.2, size = 0.5)  +
    theme_bw() +
    labs(x = "Herbivory (gC/m2/day)", y = "Piscivory (gC/m2/day)") +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))  

  rescor <- as.data.frame(rescor) %>% 
    rownames_to_column("pair")
  
  corr <- 
  ggplot(rescor) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_pointrange(aes(y = reorder(pair, desc(pair)), 
                   x = Estimate, xmin = `l-95% CI`, xmax = `u-95% CI`,
                   color = reorder(pair, (pair)))) +
    scale_color_fish_d(option = "Ostracion_whitleyi", begin = 0, end = 1) +
    scale_y_discrete(labels = rev(c("log(Nex) - log(Pex)",
                                "log(Nex) - log(Prod)",
                                "log(Nex) - log(Herb)",
                                "log(Nex) - log(Pisc)",
                                "log(Pex) - log(Prod)",
                                "log(Pex) - log(Herb)",
                                "log(Pex) - log(Pisc)",
                                "log(Prod) - log(Herb)",
                                "log(Prod) - log(Pisc)",
                                "log(Herb) - log(Pisc)"))) +
    labs(y = "", x = "Residual correlation") +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.line.x = element_line(),
          axis.ticks = element_line()) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 8))
  corr
  
  
  layout <- "
  abcd
  aefg
  hijk
  "
  plot <- 
  corr + p1 + p2 + p3 + p4 + p5 + p6 + p7 + 
    theme(axis.title.y = 
            element_text(margin = margin(r = -100, unit = "pt"))) +
    p8 +
    p9 + p10 +
    plot_layout(design = layout) + 
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave("output/plots/fig2_corplot.png", plot, width = 11, height = 7)
  ggsave("output/plots/fig2_corplot.pdf", plot, width = 8, height = 5)
  
  # supplemental graph
  loc <- as.data.frame(loc) %>% 
    rownames_to_column("pair")
  site <- as.data.frame(site) %>% 
    rownames_to_column("pair")
  
  cor_loc <- 
    ggplot(loc) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_pointrange(aes(y = reorder(pair, desc(pair)), 
                        x = Estimate, xmin = `l-95% CI`, xmax = `u-95% CI`)) +
    scale_y_discrete(labels = rev(c("log(Nex) - log(Pex)",
                                    "log(Nex) - log(Prod)",
                                    "log(Nex) - log(Herb)",
                                    "log(Nex) - log(Pisc)",
                                    "log(Pex) - log(Prod)",
                                    "log(Pex) - log(Herb)",
                                    "log(Pex) - log(Pisc)",
                                    "log(Prod) - log(Herb)",
                                    "log(Prod) - log(Pisc)",
                                    "log(Herb) - log(Pisc)"))) +
    labs(y = "", x = "Locality-level correlation") +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.line.x = element_line(),
          axis.ticks = element_line()) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"))
  cor_loc
  
  cor_site <- 
    ggplot(site) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_pointrange(aes(y = reorder(pair, desc(pair)), 
                        x = Estimate, xmin = `l-95% CI`, xmax = `u-95% CI`)) +
    scale_y_discrete(labels = rev(c("log(Nex) - log(Pex)",
                                    "log(Nex) - log(Prod)",
                                    "log(Nex) - log(Herb)",
                                    "log(Nex) - log(Pisc)",
                                    "log(Pex) - log(Prod)",
                                    "log(Pex) - log(Herb)",
                                    "log(Pex) - log(Pisc)",
                                    "log(Prod) - log(Herb)",
                                    "log(Prod) - log(Pisc)",
                                    "log(Herb) - log(Pisc)"))) +
    labs(y = "", x = "Site-level correlation") +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.line.x = element_line(),
          axis.ticks = element_line()) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"))
  cor_site
  
  plot2 <- 
  cor_loc + cor_site 
  
  ggsave("output/plots/figSx_corplot.png", plot2, width = 8, height = 4)
  ggsave("output/plots/figSx_corplot.pdf", plot2, width = 8, height = 4)
  
  
   }



