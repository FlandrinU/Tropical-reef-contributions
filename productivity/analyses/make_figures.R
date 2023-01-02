
#' Function to run all the functions needed to make figures for paper
#' 


make_figures = function(){
  

  # #Selecting covariables of interest
  # data_formodel = RLS_Management %>%
  #   dplyr::select(-c(site_code, SurveyID, No.take.multizoned,Effectiveness, NoViolence,Voice,ControlofCorruption, log10ProdB, Productivity, Biom, Prod, Country, log10Prod:SiteLongitude,Class)) %>%
  #   na.omit()
  # 
  # cor_matrix = cor(data_formodel)
  # 
  # library(corrplot)
  # col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  # corrplot(cor_matrix, method="color", col=col(200),
  #          type="upper", order="hclust", 
  #          addCoef.col = "black", # Ajout du coefficient de corrélation
  #          tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
  #          # Combiner avec le niveau de significativité
  #          sig.level = 0.01, insig = "blank", 
  #          # Cacher les coefficients de corrélation sur la diagonale
  #          diag=FALSE)
  # 
  # 
  # #Adding date
  
#Info data 
setwd(here())
RLS_info = read.table("data/RLS_transect_info.txt") %>%
  dplyr::rename(site_code = "SiteCode") %>%
  dplyr::left_join(info, by = "site_code") %>%
  dplyr::select(site_code, SurveyID)

print(paste("Minimum biomass:",min(RLS_Management$Biom)))
print(paste("Maximum biomass:", max(RLS_Management$Biom)))
print(paste("Minimum production:",min(RLS_Management$Prod)))
print(paste("Maximum production:", max(RLS_Management$Prod)))
print(paste("Minimum turnover:",min(RLS_Management$Productivity)))
print(paste("Maximum turnover:", max(RLS_Management$Productivity)))

#Specific dataframes

traits_trophic = traits %>%
  dplyr::mutate(Species = sub(" ", "_", CURRENT_SPECIES_NAME)) %>%
  dplyr::select(Species, Trophic.Level,Trophic.group,MaxLength )

RLS_prod_all_site = RLS_prod_all %>% dplyr::left_join(RLS_info, by = "SurveyID") %>% dplyr::select(c(Species, Genus,Family,Num, Size, MaxSizeTL,W, Kmax,site_code))

RLS_prod_figures = RLS_Management %>%
  left_join(RLS_prod_all_site, by = "site_code") %>%
  left_join(traits_trophic, by = "Species") %>%
  #Trophic level groups
  mutate(TrophGroup = ifelse(Trophic.Level <= 2.5, "2/2.5",
                             ifelse(Trophic.Level <= 3, "2.5/3",
                                    ifelse(Trophic.Level <= 3.5, "3/3.5",
                                           ifelse(Trophic.Level <= 4, "3.5/4",
                                                  ifelse(Trophic.Level <= 4.5, "4/4.5",
                                                         ifelse(Trophic.Level <= 5, "4.5/5"))))))) %>%
  #Regrouping diet 
  mutate(Diet = ifelse(Trophic.group %in% c("Scraping herbivore","Algal farmer","Browsing herbivore","Corallivore","Excavator"),"Herbivore",
                       ifelse(Trophic.group %in% c("Benthic invertivore"),"Invertivore",
                              ifelse(Trophic.group %in% c("Omnivore","Cleaner"),"Omnivore",
                                     ifelse(Trophic.group %in% c("Higher carnivore","Higher carnivore"),"Carnivore","Planktivore"))))) %>% 
  #Sum by Class
  group_by(Class) %>%
  dplyr::mutate(site_biom = sum(Biom)) %>%
  ungroup() %>%
  #Sum by Class and trophic group
  group_by(Class,TrophGroup) %>%
  dplyr::mutate(troph_biom = sum(Biom)) %>%
  ungroup() %>%
  #Proportional trophic groups
  dplyr::mutate(trophic_pyramid = troph_biom/site_biom,
                trophic_percentage = trophic_pyramid * 100) %>%
  #Average by class
  group_by(Class, TrophGroup) %>%
  dplyr::mutate(class_trophic = mean(trophic_percentage)) %>%
  ungroup() %>%
  group_by(site_code) %>%
  mutate(biom_by_site = sum(Biom)) %>%
  ungroup() %>%
  group_by(site_code,Diet) %>%
  mutate(diet_biom = sum(Biom),
         diet_size = mean(Size)) %>%
  ungroup() %>%
  mutate(relative_diet = (diet_biom/biom_by_site)*100)


RLS_prod_figures$Diet <- factor(RLS_prod_figures$Diet, levels=c("Herbivore","Planktivore","Omnivore","Invertivore","Carnivore"))
group.colors.diet <- c(`Herbivore` = "#359B73", `Omnivore` = "#7D22B2", `Invertivore` = "#F9AC3B", `Planktivore` = "#C1FFC1", `Carnivore` = "#AC0624")
(diet_biomass = ggplot(RLS_prod_figures %>% filter(Class != "transition") %>% dplyr::select(Class, site_code, Diet, relative_diet) %>% distinct(), aes(Diet,relative_diet,fill=Diet))+
  geom_jitter(alpha  = 0.3, size = 0.5,color='black',pch=21) +
  geom_boxplot(alpha  = 0.8) +
  scale_fill_manual(values = group.colors.diet) +
  scale_color_manual(values = group.colors.diet) +
  ylim(0,100) +
  facet_wrap(~Class,labeller=labeller(Class=c("deadzone"="Low biomass/turnover","partial"="High turnover","pristine"="High biomass")))+
  labs(x = "",
       y= "Relative biomass (%)") +
  theme_tufte()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()))

(diet_size = ggplot(RLS_prod_figures %>% filter(Class != "transition") %>% dplyr::select(Class, site_code, Diet, diet_size) %>% distinct(), aes(Diet,diet_size,fill=Diet))+
    geom_jitter(alpha  = 0.3, size = 0.5,color='black',pch=21) +
  geom_boxplot(alpha  = 0.8) +
  scale_fill_manual(values = group.colors.diet) +
  scale_color_manual(values = group.colors.diet) +
  facet_wrap(~Class,labeller=labeller(Class=c("deadzone"=" ","partial"=" ","pristine"=" ")))+
  labs(x = "",
       y= "Mean fish size (cm)") +
  theme_tufte())

ggarrange(diet_biomass, diet_size, ncol = 1,common.legend = T)
ggsave("Figures/diet_detail.pdf",width =297,height  =210, units="mm",dpi=600)

RLS_prod_figures$TrophGroup <- factor(RLS_prod_figures$TrophGroup, levels=c("2/2.5","2.5/3","3/3.5","3.5/4","4/4.5","4.5/5"))
group.colors <- c(`2/2.5` = "#144D49", `2.5/3` = "#276C69", `3/3.5` = "#73C1C4", `3.5/4` = "lightgrey", `4/4.5` = "#BF8599",`4.5/5` = "#830342")


#Trophic pyramids
RLS_prod_figures$TrophGroup <- factor(RLS_prod_figures$TrophGroup, levels=c("2/2.5","2.5/3","3/3.5","3.5/4","4/4.5","4.5/5"))
group.colors <- c(`2/2.5` = "#144D49", `2.5/3` = "#276C69", `3/3.5` = "#73C1C4", `3.5/4` = "lightgrey", `4/4.5` = "#BF8599",`4.5/5` = "#830342")
group.colors.class <- c(deadzone = "#d69d4e", partial = "#046c9a", pristine = "#C1DB60")

protection = RLS_prod_figures %>% filter(Class == "partial") %>% dplyr::select(Class, TrophGroup, class_trophic) %>% distinct() %>% mutate(class_trophic = round(class_trophic,5)) %>%
  add_row(Class = "deadzone", TrophGroup = "4.5/5",class_trophic = 0)

protection$TrophGroup <- factor(protection$TrophGroup, levels=c("2/2.5","2.5/3","3/3.5","3.5/4","4/4.5","4.5/5"))

ggplot(protection)+
  geom_bar(aes(x = TrophGroup,y = class_trophic,fill=TrophGroup),stat='identity')+
  # scale_y_continuous(breaks = brks)+
  scale_fill_manual(values = group.colors)+
  coord_flip() +
    theme_void()+
  theme(legend.position="top")

ggsave("Figures/partial_trophic.pdf")

protection_pristine = RLS_prod_figures %>% filter(Class == "pristine") %>% dplyr::mutate(Trophic.group = ifelse(Trophic.group == "higher carnivore", "Higher carnivore",Trophic.group)) %>%
  dplyr::mutate(individual_biomass = W * Num,
                total_biomass = sum(individual_biomass)) %>%
  group_by(Trophic.group) %>%
  dplyr::mutate(trophic_group_biomass = sum(individual_biomass)) %>%
  ungroup() %>% 
  dplyr::mutate(relative_biomass = (trophic_group_biomass/total_biomass)*100) %>%
  dplyr::select(Class, Trophic.group, relative_biomass) %>%
  distinct(.keep_all = T)

protection_deadzone= RLS_prod_figures %>% filter(Class == "deadzone")  %>% dplyr::mutate(Trophic.group = ifelse(Trophic.group == "higher carnivore", "Higher carnivore",Trophic.group)) %>%
  dplyr::mutate(individual_biomass = W * Num,
                total_biomass = sum(individual_biomass)) %>%
  group_by(Trophic.group) %>%
  dplyr::mutate(trophic_group_biomass = sum(individual_biomass)) %>%
  ungroup() %>% 
  dplyr::mutate(relative_biomass = (trophic_group_biomass/total_biomass)*100) %>%
  dplyr::select(Class, Trophic.group, relative_biomass) %>%
  distinct(.keep_all = T)

protection_partial = RLS_prod_figures %>% filter(Class == "partial")  %>% dplyr::mutate(Trophic.group = ifelse(Trophic.group == "higher carnivore", "Higher carnivore",Trophic.group)) %>%
  dplyr::mutate(individual_biomass = W * Num,
                total_biomass = sum(individual_biomass)) %>%
  group_by(Trophic.group) %>%
  dplyr::mutate(trophic_group_biomass = sum(individual_biomass)) %>%
  ungroup() %>% 
  dplyr::mutate(relative_biomass = (trophic_group_biomass/total_biomass)*100) %>%
  dplyr::select(Class, Trophic.group, relative_biomass) %>%
  distinct(.keep_all = T)

(pristine = ggplot(protection_pristine,aes(Trophic.group, relative_biomass))+
  geom_bar(stat="identity",aes(fill = Class))+
  scale_fill_manual(values=group.colors.class,labels=c("Low Productivity/Biomass reefs","Productive reefs","High Biomass reefs"))+
  labs(x = "",
       y="Relative biomass (%)",
       fill = "Class")+
    ylim(0,50)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=0.9))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

(deadzone = ggplot(protection_deadzone,aes(Trophic.group, relative_biomass))+
  geom_bar(stat="identity",aes(fill = Class))+
  scale_fill_manual(values=group.colors.class,labels=c("Low Productivity/Biomass reefs","Productive reefs","High Biomass reefs"))+
  labs(x = "",
       y="Relative biomass (%)",
       fill = "Class")+
  theme_minimal()+
    ylim(0,50)+
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=0.9))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

partial = ggplot(protection_partial,aes(Trophic.group, relative_biomass))+
  geom_bar(stat="identity",aes(fill = Class))+
  scale_fill_manual(values=group.colors.class,labels=c("Low Productivity/Biomass reefs","Productive reefs","High Biomass reefs"))+
  labs(x = "",
       y="Relative biomass (%)",
       fill = "Class")+
  theme_minimal()+
  ylim(0,50)+
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=0.9))

ggarrange(pristine,deadzone,partial,ncol =1, common.legend = T)
ggsave("figures/Diet_by_class.png", width = 210, height = 297, units ="mm")


# ggplot(RLS_prod_figures,aes(log10(site_biom+1),trophic_percentage,group = TrophGroup, color = TrophGroup))+
#   geom_point(shape = ".", alpha = 0.5)+
#   geom_smooth()+
#   theme_minimal()+
#   scale_color_hp_d(option = "LunaLovegood")+
#   labs(x = "Standing biomass (g/m^2) - log scale",
#        y = "Relative biomass (%)",
#        color = "Trophic position")
# 
# ggsave("figures/deadzone_trophic.pdf")

#Plot figure representing K by size for each family and for some examples 
K_plots(RLS_prod_figures)

#Plot number of predeictions for K at family, genus and species level
#Check pred_type column for precise number
  
#Plot comparison of metrics biomass/biomass production/productivity and correlation levels
plot_metrics_comparison(RLS_Management)

#Table with family, genus used in the study
used_taxonomy(RLS_prod_figures)

#World maps of our transects with management information
map_management(RLS_Management,world)

#Various representation of our biomass/productivity bi-plot
plot_classes(RLS_Management)

#Figure 2 : variable importance in our model
plot_var_imp(model_test,model_5, model_5plus, model_15,model_15plus, model_sensitivity)

#Figure 3 : representations between variables and management classes
model_prob(RLS_Management,model_test)

#Supplementary maps of sites and transects
#MTE plots
# plot_MTE(data_prepped, fam_model_K)

#K observed versus predicted
# K_observed_predicted(data_prepped)

summary(RLS_Management)
}

# 
# (p.temperature.adjusted = ggplot(MTE_fam, aes(x=logMmax, y=TempStandardLogK)) +
#   geom_point(size=1,aes(colour=Family)) +
#   scale_colour_manual(values = coli)+
#   labs(title="logMmax predictor effect", x="logMmax",y="Temperature-adjusted logK") +
#   geom_segment(aes(x=MinlogMmax, xend=MaxlogMmax, 
#                    y=Intercept + SlopeLogMmax*MinlogMmax + SlopeInvTkbFixe*MeanInvTkb, 
#                    yend=Intercept + SlopeLogMmax*MaxlogMmax + SlopeInvTkbFixe*MeanInvTkb,
#                    size="Family (random effect)",colour=Family)) +
#   geom_segment(aes(x=MinlogMmaxFixe, xend=MaxlogMmaxFixe,
#                    y=InterceptFixe + SlopelogMmaxFixe*MinlogMmaxFixe + SlopeInvTkbFixe*MeanInvTkb, 
#                    yend=InterceptFixe + SlopelogMmaxFixe*MaxlogMmaxFixe+ SlopeInvTkbFixe*MeanInvTkb,
#                    size="Fish in general (fixed effect)"),colour="black") +
#   scale_size_manual(name="Predictions", values=c("Family (random effect)"=0.7, "Fish in general (fixed effect)"=1.2))+
#   theme_bw(base_size=10) +
#   theme(
#     plot.title = element_text(color="black", size=25, face="bold.italic",hjust=0.5),
#     axis.title.x = element_text(color="black", size=14, face="bold"),
#     axis.title.y = element_text(color="black", size=14, face="bold")
#   ))
# 
# # make the mass-corrected (growth is corrected at a mean mass) plot
# 
# (p.mass.adjusted = ggplot(MTE_fam, aes(x=InvTkb, y=MassStandardLogK, colour=Family)) +
#   geom_point(size=1) +
#   scale_colour_manual(values = coli) +
#   labs(title="InvTkb predictor effect", x="InvTkb",y="Mass-adjusted logK") +
#   geom_segment(aes(x=MinInvTkb, xend=MaxInvTkb, y=Intercept + SlopeInvTkb*MinInvTkb + SlopelogMmaxFixe*MeanlogMmax, 
#                    yend=Intercept + SlopeInvTkb*MaxInvTkb + SlopelogMmaxFixe*MeanlogMmax,size="Family (random effect)")) +
#   geom_segment(aes(x=MinInvTkbFixe, xend=MaxInvTkbFixe, 
#                    y=InterceptFixe + SlopeInvTkbFixe*MinInvTkbFixe + SlopelogMmaxFixe*MeanlogMmax, 
#                    yend=InterceptFixe + SlopeInvTkbFixe*MaxInvTkbFixe+ SlopelogMmaxFixe*MeanlogMmax,
#                    size="Fish in general (fixed effect)"),colour="black") +
#   scale_size_manual(name="Predictions", values=c("Family (random effect)"=0.7, "Fish in general (fixed effect)"=1.2))+
#   theme_bw(base_size=10) +
#   theme(
#     plot.title = element_text(color="black", size=25, face="bold.italic",hjust=0.5),
#     axis.title.x = element_text(color="black", size=14, face="bold"),
#     axis.title.y = element_text(color="black", size=14, face="bold")
#   ))
