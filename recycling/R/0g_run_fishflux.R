## cleaning memory
rm(list=ls())

##### 1) Get all parameters #####

# Get a dataset with unique set of species, size and sst
load(here::here("data", "metadata_surveys.Rdata") )
load(here::here("data", "data_surveys.Rdata") )
load(here::here("data", "data_species.Rdata") )

spcombo <- metadata_surveys %>%
  dplyr::select(SurveyID, sst = SiteMeanSST) %>%
  dplyr::right_join(data_surveys, by="SurveyID") %>%
  dplyr::select(sst, species, size_class) %>%
  dplyr::mutate(sst = round(sst)) %>%
  unique()

nrow(spcombo) # 19 422
length(unique(spcombo$species)) # 1024 species present out of the 1 110 with data

# Get all parameters that are independent from sst
parameters <- read_csv(here::here("data_raw", "source", "species_parameters.csv") )
sp_par <- left_join(spcombo, parameters)


# Add metabolic parameters
 
  metpar <- read.csv(here::here("data_raw",  "source", "metpar_fam_smr.csv") )  %>%
    rename(family = Family)
  B0_mean <- mean(metpar$B0)
  B0_sd <- sd(metpar$B0)
  alpha_mean <- mean(metpar$alpha)
  alpha_sd <- sd(metpar$alpha)
  theta_mean <- mean(metpar$theta)
  
  sp_par <- sp_par %>% 
    left_join(metpar)
  
  sp_par[is.na(sp_par$alpha), "alpha"] <- alpha_mean
  sp_par[is.na(sp_par$alpha_sd), "alpha_sd"] <- alpha_sd
  sp_par[is.na(sp_par$B0), "B0"] <- B0_mean
  sp_par[is.na(sp_par$B0_sd), "B0_sd"] <- B0_sd
  sp_par[is.na(sp_par$theta), "theta"] <- theta_mean
  
   # temperature adjustment of metabolic constant 
  sp_par$B0_adj <- 
    sp_par$B0 * exp(0.59 / 8.62e-5 * (1 / (28 + 273.15) - 1 / (sp_par$sst + 273.15)))
  
  sp_par <- rename(sp_par, alpha_m = alpha,
                          f0_m = B0_adj, f0_sd = B0_sd, theta_m = theta) %>%
    select(-B0)
  

# k parameter
  
# load combined data of reef services and renato's extraction of fishbase
kmax <- read_csv( here::here("data_raw", "source", "kmax_combined.csv") ) %>%
  mutate(species = gsub(" ", "_", Species), sst = sstmean, linf_m = sizemax) 

# get fishtree
fishtree <- fishtree::fishtree_complete_phylogeny(kmax$species)
# just use one tree
set.seed(1)
tree <- fishtree[[sample(1:100, 1)]]
# correlation matrix
A <- ape::vcv(tree, cor = TRUE)

# only use names that exist in tree
kmax <- filter(kmax, species %in% colnames(A))


fit_kmax <- brm(
  log(kmax) ~ log(linf_m) + sst + (1|gr(species, cov = A)), 
  data = kmax, data2 = list(A = A),
  family = gaussian(), cores = 4
)

summary(fit_kmax)
bayes_R2(fit_kmax)
# phylogenetic signal
hyp <- "sd_species__Intercept^2 / (sd_species__Intercept^2 + sigma^2) = 0"
hypothesis(fit_kmax, hyp, class = NULL)


# extrapolate
draws <- tidybayes::spread_draws(fit_kmax, r_species[species,Intercept]) %>%
  mean_qi() 

rsp <- draws$r_species
names(rsp) <- draws$species
hist(draws$r_species)

# full tree
load( here::here("data_raw",  "source", "fishtree_glob.RData") )
# just use one tree
set.seed(2)
treeglob <- fishtree[[sample(1:100, 1)]]

rphy_pred <- picante::phyEstimate(treeglob, rsp) %>%
  rownames_to_column("species") %>%
  mutate(r_phylo = estimate) %>%
  select(species, r_phylo)
  
rphy <- data.frame(
  species = names(rsp),
  r_phylo = rsp
)

rphy <- bind_rows(rphy, rphy_pred)

get_variables(fit_kmax)

vars <- fit_kmax %>%
  spread_draws(b_Intercept, b_loglinf_m, b_sst)

kpred <- lapply(1:nrow(sp_par), function(i){
  data <-
    left_join(sp_par[i,],  rphy)
  
  linf <- unique(log(data$linf_m)) 
  sst <-  unique(data$sst)
  
  kmax <- exp(vars$b_Intercept +
    (linf * vars$b_loglinf_m) +
      (sst * vars$b_sst) +
    (data$r_phylo))
  
  data.frame(
    species = data$species,
    sst = data$sst,
    k_m = mean(kmax),
    k_sd = sd(kmax)
  )
}) %>% plyr::ldply()

sp_par <- left_join(sp_par, kpred)

# Add values for AE
ae <- read_csv(here::here("data_raw", "source", "ae_dietcat.csv")) %>%
  mutate(diet_cat = as.character(diet_cat)) %>%
  mutate(an_sd = case_when(an_sd>0.2 ~ 0.2, TRUE ~ an_sd),
         ap_sd = case_when(ap_sd>0.2 ~ 0.2, TRUE ~ ap_sd),
         ac_sd = case_when(ac_sd>0.2 ~ 0.2, TRUE ~ ac_sd)) %>%
  mutate(diet_cat = as.character(diet_cat))

# mean for missing diet category
ae_5 <- ae %>%
  select(- diet_cat) %>%
  summarize_all(mean) %>%
  mutate(diet_cat = "5")

ae <- bind_rows(ae, ae_5) 

sp_par <- left_join(
  mutate(sp_par, diet_cat = as.character(diet_cat)), ae) %>%
  mutate(v_m = sst) %>%
  unique()

nrow(sp_par)
write_csv(sp_par, file=here::here("data","parameters_sp_sst.csv")  )
write_csv(sp_par, file=here::here("recycling", "outputs", "parameters_sp_sst.csv")  )


##### 2) Run fishflux #####

sp_par <- read_csv( here::here("data", "parameters_sp_sst.csv") )

data <- sp_par 
#data <- data[1:100,]
cnpflux <- parallel::mclapply(1:nrow(data), function(x){
  print(x)
  
  dt <- data[x,] 
  par <- dt %>% select(-species, - sst, - size_class, -family, -diet_cat) %>% as.list()
  mod <- fishflux::cnp_model_mcmc(TL = dt$size_class,
                                  param = par, iter = 1000)
  
  
  extr <- fishflux::extract(mod, par = c("F0c", "F0n", "F0p", "Gc", "Gn", "Gp", "Sc", "Sn", "Sp", 
                                         "Ic", "In", "Ip", "Wc", "Wn", "Wp", "Fc", "Fn", "Fp"))
  extr <- cbind(dt[,1:5], extr) 
  lim <- fishflux::limitation(mod, plot = FALSE)
  extr$limitation <-first(lim[lim$prop_lim == max(lim$prop_lim), "nutrient"])
  
  return(extr)
}, mc.cores = 4) %>% plyr::ldply()


cnpflux <- select(cnpflux, - Qc_m, - TL)

# saving
write_csv(cnpflux, here::here("recycling","outputs", "cnpflux_sp_size_sst.csv"))
