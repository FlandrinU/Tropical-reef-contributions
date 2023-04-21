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
pkgs <- c("here", "ggplot2", "grid", "gridExtra")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(file = here::here("outputs", "NN_NS_score_wheighted_mean.Rdata"))

library(ggplot2)

## --------------- NN against NS -------------
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
             size = 3, shape = 19, 
             aes(colour= rank_u_r)) +
  scale_colour_gradient(name="up_right",
                        low = "white", high="darkgoldenrod3",
                        limits =quantile(NN_NS_with_product$rank_u_r,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(data = NN_NS_with_product[which(NN_NS_with_product$rank_u_r >=
                      quantile(NN_NS_with_product$rank_u_r, probs=c(0.95), na.rm=T)), ],
             aes(y= NS_score, x = NN_score),
             size = 3, shape = 1, stroke = 1)+
  guides(color = "none") + ggnewscale::new_scale("colour") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, up_left == 1),
             size = 3, shape = 19, 
             aes(colour= rank_u_l)) +
  scale_colour_gradient(name="up_left",
                        low = "white", high="dodgerblue3",
                        limits =quantile(NN_NS_with_product$rank_u_l,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score),
             size = 3, stroke = 1, shape = 1,
             data = NN_NS_with_product[which(NN_NS_with_product$rank_u_l >=
                                               quantile(NN_NS_with_product$rank_u_l, probs=c(0.95), na.rm=T)), ] )+
  guides(color = "none") + ggnewscale::new_scale("colour") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, down_right == 1),
             size = 3, shape = 19, 
             aes(colour= rank_d_r)) +
  scale_colour_gradient(name="down_right",
                        low = "white", high="forestgreen",
                        limits =quantile(NN_NS_with_product$rank_d_r,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score),
             size = 3, stroke = 1,shape = 1,
             data = NN_NS_with_product[which(NN_NS_with_product$rank_d_r >=
                                               quantile(NN_NS_with_product$rank_d_r, probs=c(0.95), na.rm=T)), ] )+
  guides(color = "none") + ggnewscale::new_scale("colour") +
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, down_left == 1),
             size = 3, shape = 19,
             aes(colour= rank_d_l)) +
  scale_colour_gradient(name="down_left",
                        low = "white", high="grey40",
                        limits =quantile(NN_NS_with_product$rank_d_l,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score),
             size = 3, stroke = 1, shape = 1,
             data = NN_NS_with_product[which(NN_NS_with_product$rank_d_l >=
                                               quantile(NN_NS_with_product$rank_d_l, probs=c(0.95), na.rm=T)), ] )+
  
  #Add exemples
  ggrepel::geom_label_repel(
    data = dplyr::filter(NN_NS_with_product,
                         SiteCode %in% list_sites_ggrepel,
                         up_right ==1),
    aes(label= paste(SiteEcoregion)), size=4, nudge_x = 0.2, nudge_y = 0.2)+  
  ggrepel::geom_label_repel(
      data = dplyr::filter(NN_NS_with_product,
                           SiteCode %in% list_sites_ggrepel,
                           down_right ==1),
      aes(label= paste(SiteEcoregion)), size=4, box.padding = 0.7,
      nudge_x = 0.2, nudge_y = -0.2)+
  ggrepel::geom_label_repel(
    data = dplyr::filter(NN_NS_with_product,
                         SiteCode %in% list_sites_ggrepel,
                         down_left ==1),
    aes(label= paste(SiteEcoregion)), size=4, nudge_x = -0.2, nudge_y = -0.1)+
  ggrepel::geom_label_repel(
    data = dplyr::filter(NN_NS_with_product,
                         SiteCode %in% list_sites_ggrepel,
                         up_left ==1),
    aes(label= paste(SiteEcoregion)), size=4, nudge_x = -0.2, nudge_y = 0.2)+
  

  # # see MPAs
  # scale_shape_manual(values=c(20,17,18))+
  
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
  

  xlim(c(-1.8,1.8)) +
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
        legend.position = c(0.9,0.1),
        legend.background = element_rect())+
  guides(color = "none", shape = guide_legend("Protection status")) +
  
  #add arrows in margin
  annotation_custom( 
    grob = grid::linesGrob(arrow=grid::arrow(type="closed", length=unit(3,"mm")), 
                           gp=grid::gpar(col="forestgreen", lwd=6)), 
    xmin = -1.1, xmax = -0.4, ymin =-2.17, ymax = -2.17) +
  annotation_custom( 
    grob = grid::linesGrob(arrow=grid::arrow(type="closed", length=unit(3,"mm")), 
                           gp=grid::gpar(col="dodgerblue3", lwd=6)), 
    xmin = -2.1, xmax = -2.1, ymin =-0.9, ymax = -0.2) 


NN_NS_plot

plot_with_arrows <- ggplot_gtable(ggplot_build(NN_NS_plot))
plot_with_arrows$layout$clip[gt$layout$name == "panel"] <- "off"
grid::grid.draw(plot_with_arrows)



## Stack plot mpa proportion
mpa_proportion <- data.frame(quarter=NA, protection=NA, Freq=NA)
quart = c("up_right", "up_left", "down_left", "down_right")
names = c("Top NN and NS", "NS only", "Dark NN ent NS", "NN only")
for(i in c(1:4)){
  df <- cbind(rep(names[i], 3),
              as.data.frame(table(dplyr::filter(
                NN_NS_with_product, is.na(get(quart[i]))==F)$protection)))
  colnames(df) <- colnames(mpa_proportion)
  mpa_proportion <- rbind(mpa_proportion,df)
}


dplyr::mutate(pct = Freq / sum(Freq) *100)

# Stacked
ggplot(na.omit(mpa_proportion), aes(fill=protection, y=Freq, x=quarter)) + 
  geom_bar(position="fill", stat="identity",
           aes(linewidth= 0.5, color = quarter))+
  scale_fill_grey()


