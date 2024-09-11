################################################################################
##
## From the NN and NP scores assessed in '1d_weighted_mean_NP_NN_score.R', this 
##  script reproduces the figure 2 of the paper Flandrin et al. 2023.
##
## 2b_make_fig_2.R
##
## 19/04/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
# pkgs <- c("here", "ggplot2", "grid", "gridExtra", "ggpp", "dplyr", "mdthemes",
#           "cowplot", "ggtext", "sf", "magick")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
library(ggplot2)

rm(list=ls())

# ------------loading data-------------
load(file = here::here("outputs", "NN_NP_score_wheighted_mean.Rdata"))
# load(file = here::here("outputs", "NN_NP_score_wheighted_mean_SST20.Rdata")) # make figures with tropical points only
# load(file = here::here("outputs", "NN_NP_score_wheighted_mean_RANDOM.Rdata")) # make figures with random associations between contributions


coast_pacific_centered <- sf::st_read(here::here("data", "ShapeFiles coast",
                                                 "shapefile_coast_pacific_centered.shp"))

# --------------- Figure 2a: NN against NP -------------
#### Define colors by quarter ####
NN_NP_with_product <- NN_NP_scores |>
  dplyr::bind_cols( protection = ifelse(NN_NP_scores$mpa_enforcement == "High" &
                                          stringr::str_detect(NN_NP_scores$protection_status, "No take"),
                                        "No take",ifelse(is.na(NN_NP_scores$mpa_name)==F, 
                                                         "Restricted", "Fished"))) |>
  
  dplyr::mutate(NNxNP = abs(NN_score * NP_score)) |>
  dplyr::bind_cols(rank = rank(abs(NN_NP_scores$NN_score * NN_NP_scores$NP_score))) |>
  dplyr::bind_cols(up_right = ifelse(NN_NP_scores$NN_score > 0 & NN_NP_scores$NP_score > 0,1,NA)) |>
  dplyr::bind_cols(up_left = ifelse(NN_NP_scores$NN_score < 0 & NN_NP_scores$NP_score > 0,1,NA)) |>
  dplyr::bind_cols(down_right = ifelse(NN_NP_scores$NN_score > 0 & NN_NP_scores$NP_score < 0,1,NA)) |>
  dplyr::bind_cols(down_left = ifelse(NN_NP_scores$NN_score < 0 & NN_NP_scores$NP_score < 0,1,NA))

NN_NP_with_product <- NN_NP_with_product |>
  dplyr::bind_cols(rank_u_r = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$up_right, na.last="keep")) |>
  dplyr::bind_cols(rank_u_l = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$up_left, na.last="keep")) |>
  dplyr::bind_cols(rank_d_r = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$down_right, na.last="keep")) |>
  dplyr::bind_cols(rank_d_l = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$down_left, na.last="keep"))

NN_NP_with_product$protection <- factor(NN_NP_with_product$protection , levels = c("No take", "Restricted", "Fished"))

summary(NN_NP_with_product)


#### plot NN against NP ####
NN_NP_plot <- ggplot(NN_NP_with_product, 
                     aes( y= NP_score, x = NN_score) ) +

#up right quarter
geom_point(data= dplyr::filter(NN_NP_with_product, up_right == 1),
           size = 4,  
           stroke=0,
           aes(fill= rank_u_r, shape = protection)) +
  scale_fill_gradient(name="up_right",
                      low = "grey90", high="darkgoldenrod3",
                      limits =quantile(NN_NP_with_product$rank_u_r,
                                       probs=c( 0, 1), na.rm=T),  
                      na.value=NA, breaks = seq(1,400, 10)) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NP_with_product, up_left == 1),
             size = 4, 
             stroke = 0,
             aes(fill= rank_u_l, shape = protection)) +
  scale_fill_gradient(name="up_left",
                      low = "grey90", high="dodgerblue3",
                      limits =quantile(NN_NP_with_product$rank_u_l,
                                       probs=c( 0,1), na.rm=T), na.value=NA) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NP_with_product, down_right == 1),
             size = 4, 
             stroke=0,
             aes(fill= rank_d_r, shape = protection)) +
  scale_fill_gradient(name="down_right",
                      low = "grey90", high="forestgreen",
                      limits =quantile(NN_NP_with_product$rank_d_r,
                                       probs=c( 0,1), na.rm=T), na.value=NA) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NP_with_product, down_left == 1),
             size = 4,
             stroke=0,
             aes(fill= rank_d_l, shape = protection)) +
  scale_fill_gradient(name="down_left",
                      low = "grey90", high="grey30",
                      limits =quantile(NN_NP_with_product$rank_d_l,
                                       probs=c( 0,1), na.rm=T), na.value=NA) +
  guides(fill = "none") +
  
  
  # Add outliers
  geom_point(data = NN_NP_with_product[ c(
    which(NN_NP_with_product$rank_u_r >= 
            quantile(NN_NP_with_product$rank_u_r, probs=c(0.95), na.rm=T)) ,
    which(NN_NP_with_product$rank_d_l >= 
            quantile(NN_NP_with_product$rank_d_l, probs=c(0.95), na.rm=T)) ,
    which(NN_NP_with_product$rank_d_r >=
            quantile(NN_NP_with_product$rank_d_r, probs=c(0.95), na.rm=T)) ,
    which(NN_NP_with_product$rank_u_l >=
            quantile(NN_NP_with_product$rank_u_l, probs=c(0.95), na.rm=T))
  ), ],
  aes(y= NP_score, x = NN_score, shape = protection),
  size = 4,  
  stroke = 1)+
  
  
  # point shape according to the protection status
  scale_shape_manual(values=c(24,23,21))+
  
# add 50% 'square': equation: abs(x)^p + abs(y)^p = 1 -> p=0.5
geom_function(aes(x= NN_score, y = NP_score),
              fun = function(x){(1-abs(x)^0.5)^(1/0.5) },
              xlim=c(-1,1), linetype=3, linewidth=0.5) +
  geom_function(aes(x= NN_score, y = NP_score),
                fun = function(x){-(1-abs(x)^0.5)^(1/0.5) },
                xlim=c(-1,1), linetype=3, linewidth=0.5) +
  
  #add axes
  geom_vline(xintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  geom_hline(yintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  
  xlim(c(-1.85, 1.5)) +
  ylim(c(-2, 1.7)) +
  
  labs( x=  "", y = "")+
  theme_bw(base_line_size = 0)+
  theme(
    axis.text = element_text(size=13),
    legend.position = c(0.91,0.12),
    legend.background = element_rect(colour="black", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1.3,1,1.5,1, unit='cm'))+
  guides(color = "none", shape ="none")



#### Construct gradient for legend ####
png(filename = here::here("outputs", "figures","legend_gradient_NN.png"), 
    width=40, height = 2.1, units = "cm", res = 1000)
cols = colorRampPalette(c("forestgreen", "forestgreen"))(500)
par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
plot(0, type = "n", bty = "n", xlim = c(0, length(cols)), ylim = c(-0.5, 0.5), axes = FALSE, 
     ann = FALSE)
for (i in 1:length(cols)) {
  rect(i - 1, 0.5, i, 0.5 - i /length(cols), border = NA, col = cols[i])
}
dev.off()

png(filename = here::here("outputs", "figures","legend_gradient_NP.png"), 
    width=2.8, height = 40, units = "cm", res = 500)
cols = colorRampPalette(c("dodgerblue3", "dodgerblue3"))(500)
par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
plot(0, type = "n", bty = "n", ylim = c(0, length(cols)), xlim = c(-0.5, 0.5), axes = FALSE, 
     ann = FALSE)
for (i in 1:length(cols)) {
  rect(0.5 - i /length(cols), i - 1, 0.5, i, border = NA, col = cols[i])
}
dev.off()


### Add gradient in the margins of Fig 1a & names the corners ####
grad_NN <- cowplot::ggdraw() +
  cowplot::draw_image(here::here("outputs", "figures","legend_gradient_NN.png"))
grad_NP <- cowplot::ggdraw() +
  cowplot::draw_image(here::here("outputs", "figures","legend_gradient_NP.png"))

label_corner <- data.frame(
  X = c(-Inf,-Inf,Inf,Inf),
  Y =  c(-Inf, Inf,-Inf,Inf),
  text = c("Dark spots","NP only",
           "NN only","Bright spots"),
  x_adjust = c(-0.4,-0.4,1.5,1.3),
  y_adjust = c(-2,3,-2,3),
  col = c("grey30", "dodgerblue3", "forestgreen", "darkgoldenrod3"))

fig_2a <- NN_NP_plot +
  geom_text(data=label_corner, 
            aes(x=X,y=Y, hjust= x_adjust,vjust=y_adjust,label=text, 
                fontface= "bold", color = col),
            size = 6)+
  scale_color_identity() +
  annotation_custom( ggplotGrob(grad_NN),
                     xmin= -2.05 ,#-2.2,#-2,
                     ymin = -2.8 ,#-2.5,  #-2.4,
                     xmax=1.67,  #2.2,#1,
                     ymax = -2.4 )+
  
  annotation_custom( ggplotGrob(grad_NP),
                     xmin=-2.33, #-2.7
                     xmax=-2.13, # -2.3
                     ymin = -2.2, #-2
                     ymax = 1.9)+
  labs( x=  "Nature for Nature", y = "Nature for People")+
  theme( legend.position = "none",
         axis.title.x = element_text(hjust = 0.98,
                                     vjust = -4.5,
                                     colour = "white",
                                     face = "bold",
                                     size = 16),
         axis.title.y = element_text(hjust = 0.98,
                                     vjust = 2.5,
                                     colour = "white",
                                     face = "bold",
                                     size = 16)) 


fig_2a
ggsave(plot=fig_2a, filename=here::here("outputs", "figures", "FIG2A_NNxNP_space_5percent_outliers.png"),
       height = 8.5, width = 11)


# --------------- Figure 2b: Stack plot MPA proportion -------------
mpa_proportion <- data.frame(rect = character(), tot = numeric(), 
                             quarter = character(), protection = character(),
                             Freq = integer(), pct = numeric(), 
                             stringsAsFactors = FALSE)
quart = c("up_right", "up_left", "down_left", "down_right")
names = c("4-bright", "2-NPonly", "1-dark", "3-NNonly")


for(i in c(1:4)){
  df <- cbind(c(names[i], NA,NA),
              rep(sum(NN_NP_with_product[,quart[i]], na.rm = T), 3),
              rep(names[i], 3),
              as.data.frame(table(dplyr::filter(
                NN_NP_with_product, is.na(get(quart[i]))==F)$protection))) |>
    dplyr::mutate(pct = Freq / sum(Freq) *100)
  
  colnames(df) <- colnames(mpa_proportion)
  mpa_proportion <- rbind(mpa_proportion,df)
}

#order the stack plot
mpa_proportion$protection <- factor(mpa_proportion$protection , levels = c("No take", "Restricted", "Fished"))

# Stacked
mpa_plot <- ggplot(mpa_proportion,
                   aes(x=quarter,
                       y=Freq,
                       fill= protection)) +
  geom_bar(position = position_fill(), stat = "identity", width=0.7) +
  scale_fill_grey(start=1, end=0.6) +

  geom_text(aes( label = paste0(round(pct, 0), "%")),
            stat= "identity",
            position = position_fill(vjust = 0.5),
            size=5,
            angle = 0,
            color=rep(c("black", "black", "white"), 4)) +

  
  #Add rectangles
  geom_bar(data = mpa_proportion |> dplyr::filter(!is.na(rect)),
           aes(x=rect, color=rect),
           width=0.7,
           stat = "identity",
           alpha=0, linewidth=1,
           position = position_fill())+
  scale_colour_manual(values=c("grey30",  "dodgerblue3","forestgreen", "darkgoldenrod3")) +
  # geom_rect(aes(xmin = as.numeric(factor(quarter)) - 0.35,
  #               xmax = as.numeric(factor(quarter)) + 0.35,
  #               ymin = 0,
  #               ymax = 1),
  #           fill = NA,
  #           width = 1.5)+

  scale_x_discrete(labels=c("<span style = 'color:grey30;'>**Dark spots**</span>",
                            "<span style = 'color:dodgerblue3;'>**NP only**</span>",
                            "<span style = 'color:forestgreen;'>**NN only**</span>",
                            "<span style = 'color:darkgoldenrod3;'>**Bright spots**</span>"))+
  geom_text(size = 5, fontface = "italic",
            aes(y=1.07, label = paste0("n=", tot)))+
  
  #Define theme
  labs( x=  "", y = "")+
  theme_minimal()+
  theme(axis.text.y=element_blank(), 
        axis.text.x = ggtext::element_markdown(size=13, vjust =1.2, hjust = 0.8, #vjust =0.9, hjust = 0.5,
                                               angle = 30, color="black"),
        legend.position = "none",
        legend.background = element_rect(),
        legend.text = element_text(size = 11),
        # legend.title = element_text(size = 12),
        # plot.background = element_rect(fill = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.5,0.4,-1,.2, unit = "cm"))+
  guides(color = "none", fill = guide_legend("Protection status"))

mpa_plot  



# --------------- Figure 2c Plot outliers on map (Pacific centered)-------------
#set parameters
set.seed(0621)
stroke = 0.5 
size_point = 4
width_jitter = 200000
height_jitter = 200000
qt = 0.95

#change spatial projection of points
NN_NP_with_product_spatial <- sf::st_as_sf(NN_NP_with_product,
                                           coords=c("SiteLongitude", "SiteLatitude"),
                                           crs=4326)
NN_NP_with_product_spatial <- sf::st_transform(NN_NP_with_product_spatial,
                                               crs="+proj=cea +lon_0=150 +lat_ts=30 +x_0=150 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")

NN_NP_with_product_sp <- tidyr::extract(NN_NP_with_product_spatial, geometry,
                                        into=c("SiteLongitude", "SiteLatitude"),
                                        '\\((.*),(.*)\\)', conv = T)

#order protection to keep the same classification
NN_NP_with_product_sp$protection <- factor(NN_NP_with_product_sp$protection , levels = c("No take", "Restricted", "Fished"))

#plot points on map
fig_2c <- ggplot() +
  geom_sf(data = coast_pacific_centered, color = NA, fill = "grey70") +
  
  #down left quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, stroke = stroke, 
             fill= "grey10",
             data = NN_NP_with_product_sp[
               which(NN_NP_with_product_sp$rank_d_l >=
                       quantile(NN_NP_with_product_sp$rank_d_l, probs=c(qt), na.rm=T)), ] )+
  
  
  #up left quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, 
             stroke = stroke,
             fill= "dodgerblue3",
             data = NN_NP_with_product_sp[
               which(NN_NP_with_product_sp$rank_u_l >=
                       quantile(NN_NP_with_product_sp$rank_u_l, probs=c(qt), na.rm=T)), ] )+
  
  
  #down right quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, stroke = stroke, 
             fill= "forestgreen",
             data = NN_NP_with_product_sp[
               which(NN_NP_with_product_sp$rank_d_r >=
                       quantile(NN_NP_with_product_sp$rank_d_r, probs=c(qt), na.rm=T)), ] )+
  
  
  #up right quarter
  geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
             position = position_jitter(width =width_jitter, height =height_jitter),
             size = size_point, stroke = stroke, 
             fill= "darkgoldenrod3",
             data = NN_NP_with_product_sp[
               which(NN_NP_with_product_sp$rank_u_r >=
                       quantile(NN_NP_with_product_sp$rank_u_r, probs=c(qt), na.rm=T)), ] )+
  

  
  
  # see MPAs
  scale_shape_manual(values=c(24,23,21))+
  
  coord_sf(ylim= c(-4500000,4500000),expand = FALSE) +
  theme_bw()+
  labs(x="", y= "") +
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
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.4,.6,0,0), units = , "cm")
  )

fig_2c
ggsave(plot=fig_2c, filename=here::here("outputs", "figures", "FIG2C_World_map_5percent_outliers.png"),
       height = 4, width = 12)


# ------------------------- Plot figure 2 -----------------------------
fig_2a_png <- cowplot::ggdraw() +
  cowplot::draw_image(here::here("outputs", "figures","FIG2A_NNxNP_space_5percent_outliers.png"))


## construct the legend
NN_NP_with_product$protection <- factor(NN_NP_with_product$protection , levels = c("No take", "Restricted", "Fished"))

legend_plot <- ggplot(NN_NP_with_product, 
                      aes( y= NP_score, x = NN_score) ) +
  
  geom_point(size = 5,  
             stroke=0.2,
             aes(fill= protection, shape = protection)) +
  scale_fill_grey(start=1, end=0.5) +
  scale_shape_manual(values=c(24,23,21))+
  labs( shape = "Protection status", fill = "Protection status")+
  theme_bw(base_line_size = 0)+
  theme(legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14, face="bold"),
        legend.background = element_rect(colour="black", linewidth = 0),
        legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0))+
  guides(fill = guide_legend(title.position = "top",
                             title.hjust = 0.5,
                             title.vjust = -1)) 
legend <- ggpubr::get_legend(legend_plot)


# Construct the panel with the legend
panel_leg <- gridExtra::grid.arrange(fig_2a_png, mpa_plot, fig_2c, legend,
                                     ncol=2,
                                     layout_matrix = rbind(c(1,1,1,1,1,1,2,2,2),
                                                           c(1,1,1,1,1,1,2,2,2),
                                                           c(1,1,1,1,1,1,2,2,2),
                                                           c(1,1,1,1,1,1,2,2,2),
                                                           c(1,1,1,1,1,1,4,4,4),
                                                           c(3,3,3,3,3,3,3,3,3),
                                                           c(3,3,3,3,3,3,3,3,3),
                                                           c(3,3,3,3,3,3,3,3,3)))



## Save figure 2
png(filename = here::here("outputs", "figures","FIGURE_2.png"), 
    width=33, height = 24, units = "cm", res = 1000)
grid::grid.newpage()
grid::grid.draw(panel_leg)
# Add labels 
grid::grid.text("A", x=unit(0.01, "npc"), y=unit(0.98, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
grid::grid.text("B", x=unit(0.65, "npc"), y=unit(0.98, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
grid::grid.text("C", x=unit(0.01, "npc"), y=unit(0.35, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
dev.off()

