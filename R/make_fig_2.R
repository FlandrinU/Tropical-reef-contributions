################################################################################
##
## Makes figure 2: NN vs NS scores
##
## make_fig_2.R
##
## 19/04/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "ggplot2", "grid", "gridExtra", "ggpp", "dplyr")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(file = here::here("outputs", "NN_NS_score_wheighted_mean.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

library(ggplot2)

## --------------- Figure 2a: NN against NS -------------
#### Define colors by quarter ####
NN_NS_with_product <- NN_NS_scores |>
  dplyr::bind_cols( protection = ifelse(NN_NS_scores$mpa_enforcement == "High" &
                                          stringr::str_detect(NN_NS_scores$protection_status, "No take"),
                                        "No take",
                                        ifelse(is.na(NN_NS_scores$mpa_name)==F, 
                                               "Restricted", "Fished"))) |>
  # # dplyr::bind_cols(
  #   protection = ifelse(is.na(NN_NS_scores$mpa_name)==F &
  #                      NN_NS_scores$mpa_enforcement != "Low" &
  #                      stringr::str_detect(NN_NS_scores$protection_status, "No take"),
  #   # protection = ifelse(is.na(NN_NS_scores$mpa_iucn_cat)==F &
  #   #                     ( NN_NS_scores$mpa_iucn_cat == "Ia" |
  #   #                         NN_NS_scores$mpa_iucn_cat == "II") ,
  #                      "protected","not protected")) |>
  dplyr::mutate(NNxNS = abs(NN_score * NS_score)) |>
  dplyr::bind_cols(rank = rank(abs(NN_NS_scores$NN_score * NN_NS_scores$NS_score))) |>
  dplyr::bind_cols(up_right = ifelse(NN_NS_scores$NN_score > 0 & NN_NS_scores$NS_score > 0,1,NA)) |>
  dplyr::bind_cols(up_left = ifelse(NN_NS_scores$NN_score < 0 & NN_NS_scores$NS_score > 0,1,NA)) |>
  dplyr::bind_cols(down_right = ifelse(NN_NS_scores$NN_score > 0 & NN_NS_scores$NS_score < 0,1,NA)) |>
  dplyr::bind_cols(down_left = ifelse(NN_NS_scores$NN_score < 0 & NN_NS_scores$NS_score < 0,1,NA))

NN_NS_with_product <- NN_NS_with_product |>
  dplyr::bind_cols(rank_u_r = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.last="keep")) |>
  dplyr::bind_cols(rank_u_l = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$up_left, na.last="keep")) |>
  dplyr::bind_cols(rank_d_r = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$down_right, na.last="keep")) |>
  dplyr::bind_cols(rank_d_l = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$down_left, na.last="keep"))

NN_NS_with_product$protection <- factor(NN_NS_with_product$protection , levels = c("No take", "Restricted", "Fished"))

summary(NN_NS_with_product)

list_sites_ggrepel <- c(
  #top NN
  "GOC11", "USEC5", "ETP208", "ETP229",
  #top NS
  "PAC61", "USEC24", "QLD28",
  #top NNxNS
  "SOL4", "RED5", "QLD135", "FP19",
  #dark NNxNS
  "CAN88", "LHI33", "WA55", "NI12"
)

#### plot NN against NN ####
NN_NS_plot <- ggplot(NN_NS_with_product, 
                     aes( y= NS_score, x = NN_score) ) +
  
  #up right quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, up_right == 1),
             size = 3,  
             stroke=0,
             aes(fill= rank_u_r, shape = protection)) +
  scale_fill_gradient(name="up_right",
                        low = "white", high="darkgoldenrod3",
                        limits =quantile(NN_NS_with_product$rank_u_r,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, up_left == 1),
             size = 3, 
             stroke = 0,
             aes(fill= rank_u_l, shape = protection)) +
  scale_fill_gradient(name="up_left",
                        low = "white", high="dodgerblue3",
                        limits =quantile(NN_NS_with_product$rank_u_l,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, down_right == 1),
             size = 3, 
             stroke=0,
             aes(fill= rank_d_r, shape = protection)) +
  scale_fill_gradient(name="down_right",
                        low = "white", high="forestgreen",
                        limits =quantile(NN_NS_with_product$rank_d_r,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, down_left == 1),
             size = 3,
             stroke=0,
             aes(fill= rank_d_l, shape = protection)) +
  scale_fill_gradient(name="down_left",
                        low = "white", high="grey40",
                        limits =quantile(NN_NS_with_product$rank_d_l,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  guides(fill = "none") +
  
  
  # Add outliers
  geom_point(data = NN_NS_with_product[ c(
    which(NN_NS_with_product$rank_u_r >= 
            quantile(NN_NS_with_product$rank_u_r, probs=c(0.95), na.rm=T)) ,
      which(NN_NS_with_product$rank_d_l >= 
            quantile(NN_NS_with_product$rank_d_l, probs=c(0.95), na.rm=T)) ,
    which(NN_NS_with_product$rank_d_r >=
            quantile(NN_NS_with_product$rank_d_r, probs=c(0.95), na.rm=T)) ,
    which(NN_NS_with_product$rank_u_l >=
            quantile(NN_NS_with_product$rank_u_l, probs=c(0.95), na.rm=T))
    ), ],
    aes(y= NS_score, x = NN_score, shape = protection),
    size = 3,  
    stroke = 1)+

  
    # see MPAs
  scale_shape_manual(values=c(24,23,21))+
    
  
  # # Add exemples
  # ggrepel::geom_label_repel(
  #   data = dplyr::filter(NN_NS_with_product,
  #                        SiteCode %in% list_sites_ggrepel,
  #                        up_right ==1),
  #   aes(label= paste(SiteEcoregion)), size=4, nudge_x = 0.2, nudge_y = 0.2)+
  # ggrepel::geom_label_repel(
  #   data = dplyr::filter(NN_NS_with_product,
  #                        SiteCode %in% list_sites_ggrepel,
  #                        down_right ==1),
  #   aes(label= paste(SiteEcoregion)), size=4, box.padding = 0.7,
  #   nudge_x = 0.2, nudge_y = -0.2)+
  # ggrepel::geom_label_repel(
  #   data = dplyr::filter(NN_NS_with_product,
  #                        SiteCode %in% list_sites_ggrepel,
  #                        down_left ==1),
  #   aes(label= paste(SiteEcoregion)), size=4, nudge_x = -0.2, nudge_y = -0.1)+
  # ggrepel::geom_label_repel(
  #   data = dplyr::filter(NN_NS_with_product,
  #                        SiteCode %in% list_sites_ggrepel,
  #                        up_left ==1),
  #   aes(label= paste(SiteEcoregion)), size=4, nudge_x = -0.2, nudge_y = 0.2)+
  
  
  # add 50% square: equation: abs(x)^p + abs(y)^p = 1 -> p=0.5
  geom_function(aes(x= NN_score, y = NS_score),
                fun = function(x){(1-abs(x)^0.5)^(1/0.5) },
                xlim=c(-1,1), linetype=3, linewidth=0.3) +
  geom_function(aes(x= NN_score, y = NS_score),
                fun = function(x){-(1-abs(x)^0.5)^(1/0.5) },
                xlim=c(-1,1), linetype=3, linewidth=0.3) +
  
  #add axes
  geom_vline(xintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  geom_hline(yintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  
  
  xlim(c(-2,2)) +
  ylim(c(-1.8,1.8)) +
  labs( x=  "Nature to Nature", y = "Nature to People")+
  theme_bw(base_line_size = 0)+
  theme(axis.title.x = element_text(hjust = 0.04,
                                    colour = "forestgreen", 
                                    face = "bold",
                                    size = 15),
        axis.title.y = element_text(hjust = 0.04,
                                    colour = "dodgerblue3", 
                                    face = "bold",
                                    size = 15),
        axis.text = element_text(size=13),
        legend.position = c(0.9,0.2),
        legend.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(color = "none", shape = guide_legend("Protection status")) +
  
  #add arrows in margin
  annotation_custom( 
    grob = grid::linesGrob(arrow=grid::arrow(type="closed", length=unit(3,"mm")), 
                           gp=grid::gpar(col="forestgreen", lwd=6)), 
    xmin = -0.95, xmax = -0.25, ymin =-2.18, ymax = -2.18) +
  annotation_custom( 
    grob = grid::linesGrob(arrow=grid::arrow(type="closed", length=unit(3,"mm")), 
                           gp=grid::gpar(col="dodgerblue3", lwd=6)), 
    xmin = -2.35, xmax = -2.35, ymin =-0.7, ymax = 0)



NN_NS_plot

plot_with_arrows <- ggplot_gtable(ggplot_build(NN_NS_plot))
plot_with_arrows$layout$clip[plot_with_arrows$layout$name == "panel"] <- "off"
grid::grid.draw(plot_with_arrows)




## --------------- Figure 2b: Stack plot mpa proportion -------------
mpa_proportion <- data.frame(rect= NA, tot=NA, quarter=NA, protection=NA, Freq=NA, pct=NA)
quart = c("up_right", "up_left", "down_left", "down_right")
names = c("NN +  NS +", "NN -  NS +", "NN -  NS -", "NN +  NS -")
for(i in c(1:4)){
  df <- cbind(c(names[i], NA,NA),
              rep(sum(NN_NS_with_product[,quart[i]], na.rm = T), 3),
              rep(names[i], 3),
              as.data.frame(table(dplyr::filter(
                NN_NS_with_product, is.na(get(quart[i]))==F)$protection))) |>
    dplyr::mutate(pct = Freq / sum(Freq) *100)
  
  colnames(df) <- colnames(mpa_proportion)
  mpa_proportion <- rbind(mpa_proportion,df)
}


#order the stack plot
mpa_proportion$protection <- factor(mpa_proportion$protection , levels = c("No take", "Restricted", "Fished"))
mpa_proportion$quarter <- factor(mpa_proportion$quarter , levels = c("NN -  NS -", "NN +  NS -", "NN -  NS +", "NN +  NS +"))
mpa_proportion$rect <- factor(mpa_proportion$rect , levels = c("NN -  NS -", "NN +  NS -", "NN -  NS +", "NN +  NS +"))

# Stacked
mpa_plot <- ggplot(mpa_proportion[-1,], 
                   aes(x=factor(quarter , levels = c("NN -  NS -", "NN +  NS -", "NN -  NS +", "NN +  NS +")),
                       y=Freq, 
                       fill= protection)) +
  geom_bar(position = position_fill(), stat = "identity") +
  scale_fill_grey(start=0.8, end=0.2) +
  
  geom_text(aes( label = paste0(round(pct, 0), "%")),
            stat= "identity",
            position = position_fill(vjust = 0.5),
            size=5,
            angle = 90,
            color=rep(c("black", "black", "white"), 4)) +
  
  
  #Add rectangles
  geom_bar(aes(x=factor(rect , levels = c("NN -  NS -", "NN +  NS -", "NN -  NS +", "NN +  NS +")),
               color=rect),
           stat = "identity",
           alpha=0, linewidth=1.5,
           position = position_fill())+
  scale_colour_manual(values=c("grey10", "forestgreen", "dodgerblue3", "darkgoldenrod3" )) +
  geom_text(size = 5, aes(y=1.07, label = paste0("n=", tot)))+
  
  #Define theme
  labs( x=  "", y = "")+
  theme_minimal()+
  theme(axis.text.y=element_blank(), 
        axis.text.x = element_text(size=13, face = "bold", vjust =1.2, hjust = 1,
                                   angle = 45, color="black"),
        legend.position = "none",
        legend.background = element_rect(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.5,0.1,0.5,1.5, unit = "cm"))+
  guides(color = "none", fill = guide_legend("Protection status"))

mpa_plot  



## --------------- Figure 2c: Plot outliers on map -------------
#set parameters
stroke = 0.5 
size_point = 4
width_jitter = 2
height_jitter = 2

set.seed(06)

#plot map
fig_2c <- ggplot() +
  geom_sf(data = coast, color = NA, fill = "grey70") +
  
  #down left quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, stroke = stroke, 
             fill= "grey10",
             data = NN_NS_with_product[
               which(NN_NS_with_product$rank_d_l >=
                       quantile(NN_NS_with_product$rank_d_l, probs=c(0.95), na.rm=T)), ] )+
  
  
  #up left quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, 
             stroke = stroke,
             fill= "dodgerblue3",
             data = NN_NS_with_product[
               which(NN_NS_with_product$rank_u_l >=
                       quantile(NN_NS_with_product$rank_u_l, probs=c(0.95), na.rm=T)), ] )+
  
  
  #down right quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, stroke = stroke, 
             fill= "forestgreen",
             data = NN_NS_with_product[
               which(NN_NS_with_product$rank_d_r >=
                       quantile(NN_NS_with_product$rank_d_r, probs=c(0.95), na.rm=T)), ] )+
  
  
  #up right quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, stroke = stroke, 
             fill= "darkgoldenrod3",
             data = NN_NS_with_product[
               which(NN_NS_with_product$rank_u_r >=
                       quantile(NN_NS_with_product$rank_u_r, probs=c(0.95), na.rm=T)), ] )+
  
  
  # see MPAs
  scale_shape_manual(values=c(21,24,23))+
  
  coord_sf(xlim = c(-180,180), ylim= c(-37, 32),expand = FALSE) +
  theme_bw()+
  labs(x="Longitude", y= "Latitude") +
  theme(axis.title.x = element_text(face = "bold",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    size = 15),
        axis.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey95"),
        legend.position = "none",
        plot.title = element_text(size=10, face="bold"),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
  )

fig_2c 

#### ----------------- Panel -------------------
arrange <- ggpubr::ggarrange(plot_with_arrows,
                             mpa_plot, 
                             widths = c(3,2),
                             labels = c("A", "B"),
                             font.label = list(size=17),
                             ncol = 2)
arrange

ggsave(plot=arrange, filename = here::here("outputs", "figures", "Panel_NNvsNS_and_MPA.png"),
       width = 35, height = 21 , units = "cm") #5:3 = x:y ratio





png(filename = here::here("outputs", "figures","Panel_fig_2.png"), 
    width=40, height = 27, units = "cm", res = 1000)

gridExtra::grid.arrange(
  plot_with_arrows,
  mpa_plot, 
  fig_2c,
  ncol=2,
  layout_matrix = rbind(c(1,1,1,2,2),
                        c(1,1,1,2,2),
                        c(1,1,1,2,2),
                        c(1,1,1,2,2),
                        c(3,3,3,3,3),
                        c(3,3,3,3,3),
                        c(3,3,3,3,3)))
# Add labels 
grid::grid.text("A", x=unit(0.01, "npc"), y=unit(0.98, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
grid::grid.text("B", x=unit(0.60, "npc"), y=unit(0.98, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
grid::grid.text("C", x=unit(0.01, "npc"), y=unit(0.42, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
dev.off()
