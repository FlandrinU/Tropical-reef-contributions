## cleaning memory
rm(list=ls())

library(here)
library(tidyverse)
library(ggridges)
library(patchwork)

library(ggfortify)
library(corrmorant)
library(FactoMineR)

# dataset
load(here::here("results","task3_data_sites.Rdata"))

data_plot <- task3_data_sites %>%
  mutate(type="sites") %>%
  mutate(nb_surveys = case_when ( n_surveys == 1 ~ "1",
                                  n_surveys == 2 ~ "2",
                                  n_surveys >= 3 ~ "3-6" ) 
         )

##### total biomass and density, recycling, storage

var_plot<-c("Btot_fishflux", "storage_N", "recycling_N", 
            "Ntot_fishflux", "storage_P", "recycling_P") 

for (k in 1:length(var_plot) ) {
  
  var_k<-var_plot[k]
  
  plot_k<- ggplot(data_plot, aes(x = .data[[var_k]], y = type, fill = factor(stat(quantile))) ) +
    stat_density_ridges( geom = "density_ridges_gradient", calc_ecdf = TRUE,
                         quantiles = 10, quantile_lines = TRUE) +
    scale_fill_viridis_d(name = "Deciles") +
    scale_x_log10() +
    scale_y_discrete(expand = c(0,0)) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 
  # plot_k
  if(k==1) {to_plot <- plot_k} else { to_plot <- to_plot + plot_k }
  
}# end of k

to_plot <- to_plot + 
  patchwork::plot_layout(guides="collect", ncol=3)
ggsave(here::here("plots","histograms_sites_biomabund_fluxes.png"),
       width = 15, height = 10, units = "cm")



# ratios about recycling
var_plot<-c("pexcr_recycling_N", "pexcr_recycling_P", 
            "recyc_stor_N", "recyc_stor_P")     


for (k in 1:length(var_plot) ) {
  
  var_k<-var_plot[k]
  
  plot_k<- ggplot(data_plot, aes(x = .data[[var_k]], y = type, fill = factor(stat(quantile))) ) +
    stat_density_ridges( geom = "density_ridges_gradient", calc_ecdf = TRUE,
                         quantiles = 10, quantile_lines = TRUE) +
    scale_fill_viridis_d(name = "Deciles") +
    # scale_x_continuous(trans="pseudo_log") + # scale_x_log10() +
    scale_y_discrete(expand = c(0,0)) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank()) 
  # plot_k
  if(k==1) {to_plot <- plot_k} else { to_plot <- to_plot + plot_k }
  
}# end of k

to_plot <- to_plot + 
  patchwork::plot_layout(guides="collect", ncol=2,byrow = TRUE)
ggsave(here::here("plots","histograms_sites_ratios.png"),
       width = 15, height = 10, units = "cm")




# correlation between biomass and excretion, egestion and storage
corrmorant::corrmorant( select(data_plot,
                               log10_biomass, log10_density,
                               recycling_N, recycling_P,
                               storage_N, storage_P ),
                        corr_method = "spearman", style = "blue_red", )
ggsave(here::here("plots","sites_cor_biom_recycling_storing.png"), width = 15, height = 15, units = "cm")


corrmorant::corrmorant(  select(data_plot,
                                recycling_N, recycling_P, 
                                pexcr_recycling_N, pexcr_recycling_P, 
                                recyc_stor_N, recyc_stor_P),
                         corr_method = "spearman", style = "blue_red")
ggsave(here::here("plots","sites_cor_recycling_ratios.png"), width = 20, height = 20, units = "cm")





#######
# PCA on key variables
pca_res <- data_plot %>%
  select(
         rec_N = recycling_N, rec_P=recycling_P,
         excRec_N = pexcr_recycling_N, excRec_P = pexcr_recycling_P, 
         RecSto_N = recyc_stor_N, RecSto_P=recyc_stor_P) %>%
  FactoMineR::PCA(scale.unit = TRUE )

round(pca_res$eig,2)
# => PC1and PC2 with eig >1, PC3 has eig=0.90
round(pca_res$var$contrib[,1:3],1) 

p12<- FactoMineR::plot.PCA(pca_res, choix = "var", axes = c(1,2))           
p13 <- FactoMineR::plot.PCA(pca_res, choix = "var", axes = c(3,4)) 
to_plot <- p12 + p13 + 
  patchwork::plot_layout(guides="collect", ncol=2)
ggsave(here::here("plots","sites_pca.png"),
       width = 15, height = 8, units = "cm")



# plot with correlation of recycling and storing metrics (for N) with drivers
head(task3_data_sites)

  p1a<- ggplot(task3_data_sites, aes( x = SiteMeanSST, y = recycling_N ) ) +
  geom_point(size=0.8) + geom_smooth()

  p1b <- ggplot(task3_data_sites, aes( x = SiteMeanSST, y = recyc_stor_N ) ) +
  geom_point(size=0.8) + geom_smooth()
  
  p1c <- ggplot(task3_data_sites, aes( x = SiteMeanSST, y = pexcr_recycling_N ) ) +
    geom_point(size=0.8) + geom_smooth()

  p2a <- ggplot(task3_data_sites, aes( x = size_p90, y = recycling_N ) ) +
    geom_point(size=0.8) + geom_smooth()
  
  p2b <- ggplot(task3_data_sites, aes( x = size_p90, y = recyc_stor_N ) ) +
    geom_point(size=0.8) + geom_smooth()
  
  p2c <- ggplot(task3_data_sites, aes( x = size_p90, y = pexcr_recycling_N ) ) +
    geom_point(size=0.8) + geom_smooth()
  
  
  p3a <- ggplot(task3_data_sites, aes( x = pbiom_dether, y = recycling_N ) ) +
    geom_point(size=0.8) + geom_smooth()
  
  p3b <- ggplot(task3_data_sites, aes( x = pbiom_dether, y = recyc_stor_N ) ) +
    geom_point(size=0.8) + geom_smooth()
  
  p3c <- ggplot(task3_data_sites, aes( x = pbiom_dether, y = pexcr_recycling_N ) ) +
    geom_point(size=0.8) + geom_smooth()
  

  to_plot <- p1a + p1b + p1c + 
    p2a + p2b + p2c + 
    p3a + p3b + p3c 
    patchwork::plot_layout(guides="collect", ncol=3, byrow = FALSE)
  ggsave(here::here("plots","simple_drivers.png"),
         width = 15, height = 15, units = "cm")  



  ggplot(task3_data_sites, aes( x = log10_biomass, y = recycling_P,
                                col= pbiom_plankt) ) +
    geom_point(size=0.6)
  
  
  ggplot(task3_data_sites, aes( x = pbiom_pisci, y = pexcr_recycling_P
                                 ) ) +
    geom_point(size=0.6) + geom_smooth()

  
  ggplot(task3_data_sites, aes( x = size_p90, y = recyc_stor_P) ) +
    geom_point(size=0.6) + geom_smooth()
  
  

## boxplot nb surveys vs log10 biomass
head(data_plot)
plot_k<- data_plot %>%
  ggplot( aes(x = nb_surveys, y = Btot, fill = nb_surveys ) ) +
  geom_violin( draw_quantiles = c(0.5, 0.75, 0.95),  trim=TRUE)  
# geom_boxplot()
plot_k  

head(data_plot)
plot_k<- data_plot %>%
  ggplot( aes(x = nb_surveys, y = recycling_N, fill = nb_surveys ) ) +
  #   geom_boxplot()
geom_violin( draw_quantiles = c(0.5, 0.75, 0.95),  trim=TRUE)  

plot_k 

data_plot %>% 
  filter(recycling_N >100) %>%
  group_by(nb_surveys) %>%
  summarise(n())


data_plot %>% 
  filter(Btot >50000) %>%
  group_by(nb_surveys) %>%
  summarise(n())


data_plot %>% 
  mutate(bigbiom = case_when(Btot <= 100000  ~ "inf_100kg",
                             Btot > 100000  ~ "sup_100kg")) %>%
  ggplot( aes( x = SiteLongitude, y = Btot, 
               col = nb_surveys, shape = bigbiom ) ) +
  geom_point(size=0.8)
  
  
  group_by(nb_surveys) %>%
  ggplot( aes(x = bigbiom, y = SiteLatitude, col = nb_surveys ) ) +
  geom_boxplot()
