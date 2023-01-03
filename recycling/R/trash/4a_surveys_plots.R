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
load(here::here("results","task3_data_surveys.Rdata"))

data_plot <- task3_data_surveys %>%
  mutate(type="surveys")

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
ggsave(here::here("plots","histograms_surveys_biomabund_fluxes.png"),
       width = 15, height = 10, units = "cm")



# ratios about recycling
var_plot<-c("excretion_NP", "egestion_NP",
            "pexcr_recycling_N", "pexcr_recycling_P", 
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
  patchwork::plot_layout(guides="collect", ncol=3,byrow = FALSE)
ggsave(here::here("plots","histograms_surveys_ratios.png"),
       width = 15, height = 10, units = "cm")



# correlation between biomass and excretion, egestion and storage
corrmorant::corrmorant( select(data_plot,
                               log10_biomass, log10_density,
                               excretion_N, excretion_P, 
                               egestion_N, egestion_P,
                               storage_N, storage_P ),
                        corr_method = "spearman", style = "blue_red")
ggsave(here::here("plots","surveys_cor_excr_eges_storing.png"), width = 10, height = 10)


corrmorant::corrmorant(  select(data_plot,
                                log10_biomass, 
                                recycling_N, recycling_P,
                                pexcr_recycling_N, pexcr_recycling_P, 
                                recyc_stor_N, recyc_stor_P),
                         corr_method = "spearman", style = "blue_red")
ggsave(here::here("plots","surveys_cor_recycl_storing.png"), width = 10, height = 10)





#######
# PCA on key variables
pca_res <- data_plot %>%
  select(log10_biomass, 
         recycling_N, recycling_P,
         pexcr_recycling_N, pexcr_recycling_P, 
         recyc_stor_N, recyc_stor_P) %>%
  FactoMineR::PCA(scale.unit = TRUE)


round(pca_res$eig,2)
# => PC1and PC2 with eig >1, PC3 has eig=0.92

round(pca_res$var$contrib[,1:3],1) 
# => PC1 driven by biomass and density and total recycling
# => PC2 driven by contribution of excretion to recycling
# => PC3 driven by ratio of recycling vs storing

p12<- FactoMineR::plot.PCA(pca_res, choix = "var", axes = c(1,2))           
p13 <- FactoMineR::plot.PCA(pca_res, choix = "var", axes = c(1,3)) 
to_plot <- p12 + p13 + 
  patchwork::plot_layout(guides="collect", ncol=2)
ggsave(here::here("plots","surveys_pca.png"))



task3_data_surveys %>%
  filter(recycling_N>400)

