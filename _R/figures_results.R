source("~/Code/coexist/_R/call_libraries.R")
source("~/Code/coexist/_R/math_functions.R")
library(scales)
library(ggplot2)
library(viridisLite)
library(grid)

fig_width <- 35
ngen <- 10000
nsp <- 20
n_axes <- 15

load(here::here("Array_simulations.RData"))
dir.create(here::here("outputs", glue::glue("Comparison")))

# Build dataset

Abundances_all <- data.frame(Mortality = numeric(),
                             Fecundity = numeric(),
                             Seed = numeric(),
                             Seed_r = numeric(),
                             Mod = character(),
                             Nb_obs = numeric(),
                             Species = numeric(),
                             Abundance = numeric())

Species_all <- data.frame(Mortality = numeric(),
                          Fecundity = numeric(),
                          Seed = numeric(),
                          Seed_r = numeric(),
                          Mod = character(),
                          Nb_obs = numeric(),
                          N_sp = numeric(),
                          Shannon = numeric())

for (simu in c(1:nrow(Simulations))){
  print(paste("Simu", simu))
  
  mortality <- Simulations[simu,1]
  fecundity <- Simulations[simu,2]
  seed <- Simulations[simu,3]
  seed_r <- Simulations[simu,4]
  
  for (mod in c("Perf_know", "Part_know", "Part_know_IV")){
    
    if(mod == "Perf_know"){
      
      load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
      
      Abundances_tmp <- data.frame(Mortality=rep(mortality, 20),
                                   Fecundity=rep(fecundity, 20),
                                   Mod=rep(mod, 20),
                                   Nb_obs=rep(0, 20),
                                   Seed=rep(seed, 20),
                                   Seed_r=rep(seed_r, 20),
                                   Species=c(1:20),
                                   Abundance=Abundances[ngen,])
      Abundances_all <- rbind(Abundances_all, Abundances_tmp)
      
      nsp_final <- length(which(Abundances[ngen,]!=0))
      
      df_shannon <- data.frame(Species = 1:nsp,
                               Abundance = Abundances[ngen,])%>%
        dplyr::mutate(Proportion = Abundance / sum(Abundance))%>%
        dplyr::filter(Abundance > 0)%>%
        dplyr::mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
      
      shannon <- -sum(df_shannon$prop_times_ln_prop)
      
      Species_tmp <- data.frame(Mortality=mortality,
                                Fecundity=fecundity,
                                Mod=mod,
                                Nb_obs=0,
                                Seed=seed,
                                Seed_r=seed_r,
                                N_sp=nsp_final,
                                Shannon=shannon)
      
      Species_all <- rbind(Species_all, Species_tmp)
      
    }else{
      
      for (nb_obs in c(0:15)){
        
        load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
        Abundances_tmp <- data.frame(Mortality=rep(mortality, 20),
                                     Fecundity=rep(fecundity, 20),
                                     Mod=rep(mod, 20),
                                     Nb_obs=rep(nb_obs, 20),
                                     Seed=rep(seed, 20),
                                     Seed_r=rep(seed_r, 20),
                                     Species=c(1:20),
                                     Abundance=Abundances[ngen,])
        Abundances_all <- rbind(Abundances_all, Abundances_tmp)
        
        nsp_final <- length(which(Abundances[ngen,]!=0))
        
        df_shannon <- data.frame(Species = 1:nsp,
                                 Abundance = Abundances[ngen,])%>%
          dplyr::mutate(Proportion = Abundance / sum(Abundance))%>%
          dplyr::filter(Abundance > 0)%>%
          dplyr::mutate(ln_prop = log(Proportion), prop_times_ln_prop = ln_prop*Proportion)
        
        shannon <- -sum(df_shannon$prop_times_ln_prop)
        
        Species_tmp <- data.frame(Mortality=mortality,
                                  Fecundity=fecundity,
                                  Mod=mod,
                                  Nb_obs=nb_obs,
                                  Seed=seed,
                                  Seed_r=seed_r,
                                  N_sp=nsp_final,
                                  Shannon=shannon)
        
        Species_all <- rbind(Species_all, Species_tmp)
      }
    }
  }
}

Species_all$Hill_Shannon <- exp(Species_all$Shannon)

save(Abundances_all, file = here::here("outputs", "Comparison", "Abundances_all.RData"))
save(Species_all, file = here::here("outputs", "Comparison", "Species_all.RData"))

Percentage_similarity <- data.frame(
  Mortality = numeric(),
  Fecundity = numeric(),
  Seed = numeric(),
  Seed_r = numeric(),
  Mod_comp = character(),
  Nb_obs = numeric(),
  PS = numeric())

for (simu in c(1:nrow(Simulations))){
  
  print(paste("Simu", simu))
  
  mortality <- Simulations[simu,1]
  fecundity <- Simulations[simu,2]
  seed <- Simulations[simu,3]
  seed_r <- Simulations[simu,4]
  
  load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
  Abundances_perf <- Abundances[ngen,]
  
  Percentage_similarity <- rbind(Percentage_similarity,
                                 data.frame(
                                   Mortality = mortality,
                                   Fecundity = fecundity,
                                   Seed = seed,
                                   Seed_r = seed_r,
                                   Mod_comp = "Perf_know",
                                   Nb_obs = NA,
                                   PS = 1
                                 ))
  
  for(nb_obs in 0:15){
    
    load(here::here("outputs", "outputs_cluster", glue::glue("Part_know_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
    Abundances_part <- Abundances[ngen,]
    
    load(here::here("outputs", "outputs_cluster", glue::glue("Part_know_IV_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
    Abundances_part_IV <- Abundances[ngen,]
    
    Percentage_similarity <- rbind(Percentage_similarity,
                                   data.frame(
                                     Mortality = rep(mortality,2),
                                     Fecundity = rep(fecundity,2),
                                     Seed = rep(seed,2),
                                     Seed_r = rep(seed_r,2),
                                     Mod_comp = c("Part_know", "Part_know_IV"),
                                     Nb_obs = rep(nb_obs,2),
                                     PS = c(percentage_similarity(Abundances_perf, Abundances_part), percentage_similarity(Abundances_perf, Abundances_part_IV))
                                   ))
  }
}

save(Percentage_similarity, file = here::here("outputs", "Comparison", "Percentage_similarity.RData"))


#### Performance on sites where the species is present in the final community ####

Perf_on_sites <- data.frame(
  Mortality = numeric(),
  Fecundity = numeric(),
  Seed = numeric(),
  Seed_r = numeric(),
  Mod = character(),
  Nb_obs = numeric(),
  Perf = numeric()
)

for (simu in c(1:nrow(Simulations))){
  print(paste("Simu", simu))
  
  mortality <- Simulations[simu,1]
  fecundity <- Simulations[simu,2]
  seed <- Simulations[simu,3]
  seed_r <- Simulations[simu,4]
  
  #Load the species performance once per configuration
  if(seed_r==1){
    load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_perf_E_Sp.RData")))
  }
  
  for (mod in c("Perf_know", "Part_know", "Part_know_IV")){
    
    if(mod == "Perf_know"){
      
      load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
      
      perf_present <- rep(NA, length(community_end))
      w0 <- (as.vector(t(community_end))==0)
      perf_present[w0] <- NA
      perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
      
      Perf_on_sites_temp <- data.frame(
        Mortality = mortality,
        Fecundity = fecundity,
        Seed = seed,
        Seed_r = seed_r,
        Mod = mod,
        Nb_obs = NA,
        Perf = round(mean(perf_present, na.rm = TRUE), digits = 2)
      )
      
      Perf_on_sites <- rbind(Perf_on_sites, Perf_on_sites_temp)
      
    }else{
      
      for (nb_obs in c(0:15)){
        
        load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
        
        perf_present <- rep(NA, length(community_end))
        w0 <- (as.vector(t(community_end))==0)
        perf_present[w0] <- NA
        perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
        
        Perf_on_sites_temp <- data.frame(
          Mortality = mortality,
          Fecundity = fecundity,
          Seed = seed,
          Seed_r = seed_r,
          Mod = mod,
          Nb_obs = nb_obs,
          Perf = round(mean(perf_present, na.rm = TRUE), digits = 2)
        )
        
        Perf_on_sites <- rbind(Perf_on_sites, Perf_on_sites_temp)
        
      }#for n_obs
    }#else Part_know
  }#for mod
}#for simu

save(Perf_on_sites, file = here::here("outputs", "Comparison", "Perf_on_sites.RData"))

# Retrieve R2 of statistical model #

R2_df <- data.frame(Mortality=character(),
                    Fecundity=character(),
                    Seed=integer(),
                    Nb_obs=integer(),
                    R2=numeric())

IV_all <- data.frame(Mortality=character(),
                     Fecundity=character(),
                     Seed=integer(),
                     Nb_obs=integer(),
                     IV=numeric())

for (simu in c(1:nrow(Simulations[Simulations[,4]==1,]))){
  
  print(paste("Simu", simu))
  
  mortality <- Simulations[Simulations[,4]==1,][simu,1]
  fecundity <- Simulations[Simulations[,4]==1,][simu,2]
  seed <- Simulations[Simulations[,4]==1,][simu,3]
  
  load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_perf_E_Sp.RData")))
  load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_env.RData")))
  load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_sites.RData")))
  load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_niche_optimum.RData")))
  
  nsp <- ncol(perf_E_Sp)
  n_axes <- length(env)
  
  for(n_observed_axes in c(0:n_axes)){
    
    if(n_observed_axes>0){
      Obs_env <- sites[,1:n_observed_axes]
    }
    
    # Data-set
    df <- data.frame(perf_E_Sp)
    names(df) <- c(sprintf("Sp %02d", 1:nsp))
    
    df_perf <- data.frame(matrix(nrow=nrow(df), ncol=ncol(df)+2*n_axes))
    df_perf[,1:ncol(df)] <- df
    for(k in 1:n_axes){
      df_perf[,ncol(df)+k] <- raster::values(raster::raster(env[[k]]))
      df_perf[,ncol(df)+(k+n_axes)] <- (df_perf[,ncol(df)+k])^2
    }
    
    colnames(df_perf) <- c(sprintf("Sp %02d", 1:nsp), sprintf("Env_%d", 1:n_axes), sprintf("Env_%d_sq", 1:n_axes))
    
    df_perf <- df_perf %>%
      tidyr::pivot_longer(cols=c("Sp 01":glue::glue("Sp {nsp}")), names_to="Species", values_to="Perf")
    
    
    if(n_observed_axes>0){
      
      formula <- as.formula(paste0("Perf~-1+Species+Species:",
                                   paste0(colnames(df_perf)[1:n_observed_axes], collapse= "+Species:"),
                                   "+Species:",
                                   paste0(colnames(df_perf)[(n_axes+1):(n_axes+n_observed_axes)], collapse= "+Species:")))
    }
    if(n_observed_axes==0){
      formula <- as.formula(paste0("Perf~-1+Species"))
    }
    
    lm_fit <- lm(formula, data=df_perf)
    
    R2_df_temp <- data.frame(Mortality=mortality,
                             Fecundity=fecundity,
                             Seed=seed,
                             Nb_obs=n_observed_axes,
                             R2=summary(lm_fit)$r.squared)
    R2_df <- rbind(R2_df, R2_df_temp)
    
    load(here::here("outputs", "outputs_cluster", glue::glue("V_intra_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_1.RData")))
    IV_all_temp <- data.frame(Mortality=rep(mortality, nsp),
                              Fecundity=rep(fecundity, nsp),
                              Seed=rep(seed, nsp),
                              Nb_obs=rep(n_observed_axes, nsp),
                              IV=V_intra$V)
    IV_all <- rbind(IV_all, IV_all_temp)
    
  }
}

save(R2_df, file = here::here("outputs", "Comparison", "R2_df.RData"))
save(IV_all, file = here::here("outputs", "Comparison", "IV_all.RData"))

## PLOTS ##

source("~/Code/coexist/_R/call_libraries.R")
source("~/Code/coexist/_R/math_functions.R")

fig_width <- 35

load(here::here("Array_simulations.RData"))
ngen = 10000
nsp = 20
n_axes = 15

load(here::here("outputs", "Comparison", "Species_all.RData"))
load(here::here("outputs", "Comparison", "Percentage_similarity.RData"))
load(here::here("outputs", "Comparison", "Perf_on_sites.RData"))
load(here::here("outputs", "Comparison", "R2_df.RData"))
load(here::here("outputs", "Comparison", "IV_all.RData"))

Summary_level_explanation_axes_nb <- R2_df%>%
  dplyr::group_by(Mortality, Fecundity, Nb_obs)%>%
  dplyr::mutate(Mean_explanation=mean(R2), Sd=sd(R2))%>%
  dplyr::slice(1)%>%
  dplyr::ungroup()%>%
  dplyr::select(-Seed, -R2)

#One plot per Mortality * Fecundity option
for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")) {
    
    # Level of IV and of explained variation #
    
    data_figure <- IV_all[which(IV_all$Mortality==mortality&IV_all$Fecundity==fecundity),]
    
    data_figure_2 <- Summary_level_explanation_axes_nb[which(Summary_level_explanation_axes_nb$Mortality==mortality&Summary_level_explanation_axes_nb$Fecundity==fecundity),]
    
    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=IV))+
      ggplot2::geom_ribbon(data=data_figure_2, ggplot2::aes(x=as.factor(Nb_obs), y=Mean_explanation, ymin=Mean_explanation-Sd, ymax=Mean_explanation+Sd, group=1), colour="deeppink3", fill="hotpink3", alpha=0.3)+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6)+
      ggplot2::geom_boxplot(alpha=0.6)+
      ggplot2::geom_point(data=data_figure_2, ggplot2::aes(x=as.factor(Nb_obs), y=Mean_explanation), colour="deeppink3")+
      ggplot2::geom_line(data=data_figure_2, ggplot2::aes(x=as.factor(Nb_obs), y=Mean_explanation, group=1), colour="deeppink3")+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = "Number of observed dimensions",
                    y = "Observed uIV")+
      ggplot2::scale_x_discrete(labels=c(0:n_axes))+
      ggplot2::theme(text = ggplot2::element_text(size = 20), legend.position = "none")+
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . * 1 / 1 , name = "Proportion of variance \n explained by the dimensions"))
    
    ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("IV_nb_axes_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    # Species richness and Hill-Shannon diversity index #
    
    data_figure <- Species_all[which(Species_all$Mortality==mortality&Species_all$Fecundity==fecundity),]
    
    data_figure <- data_figure[data_figure$Mod!="Perf_know",]
    
    #Compute delta between with and without uIV
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta_SR = N_sp - dplyr::lag(N_sp),
                    Delta_Hill_Shannon = Hill_Shannon - dplyr::lag(Hill_Shannon),
                    ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta_SR)==FALSE)%>%
      dplyr::ungroup()%>%
      dplyr::select(Seed, Seed_r, ID_delta, Delta_SR, Delta_Hill_Shannon)
    
    p5 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta_SR))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Species richness")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p5, filename=here::here("outputs", "Comparison", glue::glue("Delta_species_richness_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    p6 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta_Hill_Shannon))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Hill-Shannon index")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p6, filename=here::here("outputs", "Comparison", glue::glue("Delta_Hill_Shannon_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    #without deltas
    data_figure <- Species_all[which(Species_all$Mortality==mortality&Species_all$Fecundity==fecundity),]
    data_figure <- data_figure[data_figure$Mod!="Part_know",]
    data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
    data_figure$Nb_obs <- factor(data_figure$Nb_obs,
                                 levels = c(as.character(c(0:n_axes)), "PK"))
    
    data_figure$Boxplot_colour <- rep(0, nrow(data_figure))
    data_figure[which(data_figure$Mod=="Perf_know"),]$Boxplot_colour <- 1
    color_text <- c(rep("grey20", length(0:n_axes)), "darkred")
    
    p1 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=N_sp))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
      ggplot2::scale_colour_manual(values=c("black", "darkred"))+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = "Species richness with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p1, filename=here::here("outputs", "Comparison", glue::glue("Species_richness_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    p2 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Hill_Shannon))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
      ggplot2::scale_colour_manual(values=c("black", "darkred"))+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = "Hill-Shannon index with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p2, filename=here::here("outputs", "Comparison", glue::glue("Hill_Shannon_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    # Percentage similarity #
    
    data_figure <- Percentage_similarity[which(Percentage_similarity$Mortality==mortality&Percentage_similarity$Fecundity==fecundity),]
    
    #Compute delta between with and without uIV
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta= PS - lag(PS), ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta)==FALSE)%>%
      dplyr::select(Nb_obs, Seed, Seed_r, ID_delta, Delta)
    
    p7 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Similarity with perfect knowledge")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p7, filename=here::here("outputs", "Comparison", glue::glue("Delta_percentage_similarity_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    #without deltas
    data_figure <- Percentage_similarity[which(Percentage_similarity$Mortality==mortality&Percentage_similarity$Fecundity==fecundity),]
    data_figure <- data_figure[data_figure$Mod!="Part_know",]
    
    # data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
    # data_figure$Nb_obs <- factor(data_figure$Nb_obs,
    #                              levels = c(as.character(c(0:n_axes)), "PK"))
    
    color_text <- c(rep("grey20", length(0:n_axes)), "darkred")
    
    data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- 16
    data_figure$Nb_obs <- as.integer(data_figure$Nb_obs)
    data_1 <- data_figure[data_figure$Mod!="Perf_know",]
    data_2 <- data_figure[data_figure$Mod=="Perf_know",]
    
    p3 <- ggplot2::ggplot(data=data_1, ggplot2::aes(x=Nb_obs, y=PS))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::geom_boxplot(alpha=0.6, colour="black", ggplot2::aes(group=Nb_obs))+
      ggplot2::geom_point(data=data_2, ggplot2::aes(x=Nb_obs, y=PS), colour="darkred")+
      ggplot2::scale_x_continuous(breaks=c(0:(n_axes+1)), labels=c(0:n_axes, "PK"))+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = "Similarity between PK and IK with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p3, filename=here::here("outputs", "Comparison", glue::glue("Percentage_similarity_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    # Mean performance (from the perfect knowledge model) of the final species community #
    
    data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality==mortality&Perf_on_sites$Fecundity==fecundity),]
    
    data_figure <- data_figure[data_figure$Mod!="Perf_know",]
    
    #Compute delta between with and without uIV
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta = Perf - dplyr::lag(Perf) , ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta)==FALSE)%>%
      dplyr::ungroup()%>%
      dplyr::select(Seed, Seed_r, ID_delta, Delta)
    
    p8 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Mean performance")))+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14, colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p8, filename=here::here("outputs", "Comparison", glue::glue("Delta_performance_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality==mortality&Perf_on_sites$Fecundity==fecundity),]
    data_figure <- data_figure[data_figure$Mod!="Part_know",]
    data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
    data_figure$Nb_obs <- factor(data_figure$Nb_obs,
                                 levels = c(as.character(c(0:n_axes)), "PK"))
    data_figure$Boxplot_colour <- rep(0, nrow(data_figure))
    data_figure[which(data_figure$Mod=="Perf_know"),]$Boxplot_colour <- 1
    color_text <- c(rep("grey20", length(0:n_axes)), "darkred")
    
    #without deltas
    p4 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Perf))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
      ggplot2::scale_colour_manual(values=c("black", "darkred"))+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = "Mean performance with uIV")+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=14),
                     axis.text.x = ggplot2::element_text(colour = color_text),
                     axis.text.y = ggplot2::element_text(colour = "grey20"),
                     legend.position = "none")
    
    ggplot2::ggsave(p4, filename=here::here("outputs", "Comparison", glue::glue("Performance_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
    
    
    data_level_explanation <- Summary_level_explanation_axes_nb[Summary_level_explanation_axes_nb$Mortality==mortality&Summary_level_explanation_axes_nb$Fecundity==fecundity,]
    
    fun <- function(x) ifelse(x == 0, "0", sub("^0+", "", x))
    
    labels <- fun(round(data_level_explanation$Mean_explanation, digits=2))
    
    plot_level_explanation <- ggplot2::ggplot(data_level_explanation, ggplot2::aes(Nb_obs, Mean_explanation))+
      ggplot2::geom_blank()+
      ggplot2::theme_classic()+
      ggplot2::scale_x_continuous(name="Level of explained variance", breaks=data_level_explanation$Nb_obs, labels=labels)+
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size=13),
                     axis.line.x = ggplot2::element_line(arrow = grid::arrow(length=grid::unit(0.15, "inches")), colour="deeppink3"),
                     axis.ticks.x=ggplot2::element_line(colour="deeppink3"),
                     axis.text.x=ggplot2::element_text(colour="deeppink3"),
                     axis.line.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     panel.grid.minor.y=ggplot2::element_blank(),
                     panel.grid.major.y=ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(l = 50,
                                                   r = 0,
                                                   t = 0,
                                                   b = 0.5))
    
    arrange_SR_Hill_Shannon <- ggpubr::ggarrange(p1, p2, p5, p6,
                                            nrow=2, ncol=2, align = "v",
                                            labels=c("A", "B", "C", "D"))
    arrange_PS_Perf <- ggpubr::ggarrange(p3, p4, p7, p8,
                                         nrow=2, ncol=2, align = "v",
                                         labels=c("A", "B", "C", "D"))
    arrange_arrow <- ggpubr::ggarrange(plot_level_explanation, plot_level_explanation,
                                       nrow=1, ncol=2, align="v")
    arrange_results_1 <- ggpubr::ggarrange(arrange_SR_Hill_Shannon, arrange_arrow, nrow=2, ncol=1, heights=c(15, 1))
    arrange_results_2 <- ggpubr::ggarrange(arrange_PS_Perf, arrange_arrow, nrow=2, ncol=1, heights=c(15, 1))
    ggplot2::ggsave(arrange_results_1, filename=here::here("outputs", "Comparison", glue::glue("Results_1_{mortality}_{fecundity}.png")),
                    width=40, height=30, units="cm", dpi=300)
    ggplot2::ggsave(arrange_results_2, filename=here::here("outputs", "Comparison", glue::glue("Results_2_{mortality}_{fecundity}.png")),
                    width=40, height=30, units="cm", dpi=300)
    
    arrange_a <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
                                   nrow=2, ncol=4, align = "v")
    arrange_b <- ggpubr::ggarrange(plot_level_explanation, plot_level_explanation, plot_level_explanation, plot_level_explanation,
                                   nrow=1, ncol=4, align="v")
    final_plot <- ggpubr::ggarrange(arrange_a, arrange_b, nrow=2, ncol=1, heights=c(15, 1))

    ggplot2::ggsave(final_plot, filename=here::here("outputs", "Comparison", glue::glue("Results_all_{mortality}_{fecundity}.png")),
                    width=50, height=50/2, units="cm", dpi=300)
    
  }
}

### Appendix ###

## Understanding why delta performance is higher than 0 with the proportional mortality ##

Comparison_suboptimals <- data.frame(
  Mortality=character(),
  Fecundity=character(),
  Mod=character(),
  Nb_obs=integer(),
  Seed=integer(),
  Seed_r=integer(),
  Abund_suboptimals=numeric(),
  perf_suboptimals=numeric(),
  perf_not_suboptimals=numeric(),
  perf_tot=numeric())

for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")){
    
    for(seed in 1:10){
      
      #same perf_E_Sp for all repetitions of the same configuration
      load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_1_perf_E_Sp.RData")))
      winner <- apply(X=perf_E_Sp, MARGIN=1, FUN=which.max)
      
      for(seed_r in 1:10){
        mod <- "Perf_know"
        load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
        
        perf_present <- rep(NA, length(community_end))
        w0 <- (as.vector(t(community_end))==0)
        perf_present[w0] <- NA
        perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
        
        Comparison_suboptimals_temp <- data.frame(
          Mortality=mortality,
          Fecundity=fecundity,
          Mod=mod,
          Nb_obs=NA,
          Seed=seed,
          Seed_r=seed_r,
          Abund_suboptimals=length(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
          perf_suboptimals=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
          perf_not_suboptimals=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)]),
          perf_tot=mean(perf_present, na.rm = TRUE))
        
        Comparison_suboptimals <- rbind(Comparison_suboptimals, Comparison_suboptimals_temp)
        
        for (mod in c("Part_know", "Part_know_IV")){
          for (nb_obs in 1:n_axes){
            
            load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
            
            
            perf_present <- rep(NA, length(community_end))
            w0 <- (as.vector(t(community_end))==0)
            perf_present[w0] <- NA
            perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community_end))[!w0]])
            
            Comparison_suboptimals_temp <- data.frame(
              Mortality=mortality,
              Fecundity=fecundity,
              Mod=mod,
              Nb_obs=nb_obs,
              Seed=seed,
              Seed_r=seed_r,
              Abund_suboptimals=length(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
              perf_suboptimals=mean(perf_present[which((winner==as.vector(t(community_end)))==FALSE)]),
              perf_not_suboptimals=mean(perf_present[which((winner==as.vector(t(community_end)))==TRUE)]),
              perf_tot=mean(perf_present, na.rm = TRUE))
            
            Comparison_suboptimals <- rbind(Comparison_suboptimals, Comparison_suboptimals_temp)
          }
        }
      }
    }
  }
}

save(Comparison_suboptimals, file = here::here("outputs", "Comparison", "Comparison_suboptimals.RData"))

load(here::here("outputs", "Comparison", "Comparison_suboptimals.RData"))


for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")){
    
    data_figure <- Comparison_suboptimals[which(Comparison_suboptimals$Mortality==mortality&Comparison_suboptimals$Fecundity==fecundity),]
    
    data_figure <- data_figure[data_figure$Mod!="Perf_know",]
    
    #Compute delta between with and without
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta = Abund_suboptimals - dplyr::lag(Abund_suboptimals) , ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta)==FALSE)%>%
      dplyr::ungroup()%>%
      dplyr::select(Seed, Seed_r, ID_delta, Delta)
    
    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " Sub-optimal species abundance")))+
      ggplot2::ggtitle(paste0("Mortality ", mortality, ", Fecundity ", fecundity))+
      ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                     axis.text = ggplot2::element_text(size=7),
                     legend.position = "none")
    
    ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_abundance_suboptimals_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}


mortality <- c("prop", "stocha")
fecundity <- "abund"

data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality%in%mortality&Perf_on_sites$Fecundity==fecundity),]

data_figure <- data_figure[data_figure$Mod!="Perf_know",]

#Compute delta between with and without
data_figure <- data_figure%>%
  dplyr::group_by(Mortality, Seed, Seed_r, Nb_obs)%>%
  dplyr::mutate(Delta = Perf - dplyr::lag(Perf) , ID_delta=Nb_obs)%>%
  dplyr::filter(is.na(Delta)==FALSE)%>%
  dplyr::ungroup()%>%
  dplyr::select(Mortality, Seed, Seed_r, ID_delta, Delta)

p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
  ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                    shape=as.factor(Mortality),
                                    group=Mortality),
                       alpha=0.6,
                       position=ggplot2::position_jitterdodge(jitter.width=0.5),
                       show.legend=F)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
  ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"), labels=c("Deterministic", "Stochastic"))+
  ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                y = expression(paste(Delta, " Mean performance")))+
  ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size=7))

ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_performance_comparison_mortality.png")),
                width=fig_width, height=fig_width/2, units="cm", dpi=300)


mortality <- c("prop", "stocha")
fecundity <- "abund"

for (mod in c("Part_know", "Part_know_IV")){
  
  data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality%in%mortality&Perf_on_sites$Fecundity==fecundity),]
  
  data_figure <- data_figure[data_figure$Mod!="Perf_know",]
  
  data_figure <- data_figure[data_figure$Mod==mod,]
  
  p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Perf))+
    ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                      shape=as.factor(Mortality),
                                      group=Mortality),
                         alpha=0.6,
                         position=ggplot2::position_jitterdodge(jitter.width=0.5),
                         show.legend=F)+
    ggplot2::scale_colour_viridis_d()+
    ggnewscale::new_scale("colour")+
    ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
    ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
    ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                  y = "Mean performance")+
    ggplot2::ggtitle(paste0(mod))+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                   axis.text = ggplot2::element_text(size=7))
  
  ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Performance_comparison_mortality_{mod}.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

# Compare delta abundance sub-optimal species with deterministic and stochastic mortality #

mortality <- c("prop", "stocha")
fecundity <- "abund"

data_figure <- Comparison_suboptimals[which(Comparison_suboptimals$Mortality%in%mortality&Comparison_suboptimals$Fecundity==fecundity),]

data_figure <- data_figure[data_figure$Mod!="Perf_know",]

#Compute delta between with and without
data_figure <- data_figure%>%
  dplyr::group_by(Mortality, Seed, Seed_r, Nb_obs)%>%
  dplyr::mutate(Delta = Abund_suboptimals - dplyr::lag(Abund_suboptimals) , ID_delta=Nb_obs)%>%
  dplyr::filter(is.na(Delta)==FALSE)%>%
  dplyr::ungroup()%>%
  dplyr::select(Mortality, Seed, Seed_r, ID_delta, Delta)

p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
  ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                    shape=as.factor(Mortality),
                                    group=Mortality),
                       alpha=0.6,
                       position=ggplot2::position_jitterdodge(jitter.width=0.5),
                       show.legend=F)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
  ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
  ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                y = expression(paste(Delta, " Sub-optimal species abundance")))+
  ggplot2::theme(text = ggplot2::element_text(size = 20),
                 axis.text = ggplot2::element_text(size=16),
                 legend.position = "none")

ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_abundance_suboptimals_comparison_mortality.png")),
                width=fig_width, height=fig_width/2, units="cm", dpi=300)

# Sub-optimal species abundance with deterministic and stochastic mortality #

mortality <- c("prop", "stocha")
fecundity <- "abund"

for (mod in c("Perf_know", "Part_know", "Part_know_IV")){
  
  data_figure <- Comparison_suboptimals[which(Comparison_suboptimals$Mortality%in%mortality&Comparison_suboptimals$Fecundity==fecundity),]
  
  data_figure <- data_figure[data_figure$Mod==mod,]
  
  if(mod=="Perf_know"){
    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Abund_suboptimals))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                        shape=as.factor(Mortality),
                                        group=Mortality),
                           alpha=0.6,
                           position=ggplot2::position_jitterdodge(jitter.width=0.5),
                           show.legend=F)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
      ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
      ggplot2::labs(x="", y = "Sub-optimal species abundance")+
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     axis.text = ggplot2::element_text(size=16),
                     axis.text.x=ggplot2::element_blank(), #remove x axis labels
                     axis.ticks.x=ggplot2::element_blank(),
                     legend.position = "none")
  }else{
    
    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Abund_suboptimals))+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed),
                                        shape=as.factor(Mortality),
                                        group=Mortality),
                           alpha=0.6,
                           position=ggplot2::position_jitterdodge(jitter.width=0.5),
                           show.legend=F)+
      ggplot2::scale_colour_viridis_d()+
      ggnewscale::new_scale("colour")+
      ggplot2::geom_boxplot(ggplot2::aes(colour=Mortality), alpha=0.6)+
      ggplot2::scale_colour_manual(values=c("#808080","black", "#4C5866"))+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = "Sub-optimal species abundance")+
      ggplot2::theme(text = ggplot2::element_text(size = 20),
                     axis.text = ggplot2::element_text(size=16),
                     legend.position = "none")
  }
  
  ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Abund_suboptimals_comparison_mortality_{mod}.png")),
                  width=fig_width, height=fig_width/2, units="cm", dpi=300)
}

# Delta mean performance vs. delta abundance of sub-optimal species #

fecundity <- "abund"
for (mortality in c("prop", "stocha")){
  
  Perf_vs_suboptimal_abundance <- Comparison_suboptimals[which(Comparison_suboptimals$Mortality==mortality&Comparison_suboptimals$Fecundity==fecundity),]
  Perf_vs_suboptimal_abundance <- Perf_vs_suboptimal_abundance[Perf_vs_suboptimal_abundance$Mod!="Perf_know",]
  
  Perf_vs_suboptimal_abundance <- Perf_vs_suboptimal_abundance%>%
    dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
    dplyr::mutate(Delta_perf = perf_tot - dplyr::lag(perf_tot) , ID_delta=Nb_obs)%>%
    dplyr::mutate(Delta_abund_suboptimal = Abund_suboptimals - dplyr::lag(Abund_suboptimals) , ID_delta=Nb_obs)%>%
    dplyr::filter(is.na(Delta_perf)==FALSE, is.na(Delta_abund_suboptimal)==FALSE)%>%
    dplyr::ungroup()%>%
    dplyr::select(Seed, Seed_r, ID_delta, Delta_perf, Delta_abund_suboptimal)
  
  p <- ggplot2::ggplot(data=Perf_vs_suboptimal_abundance, ggplot2::aes(x=Delta_abund_suboptimal, y=Delta_perf))+
    ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
    ggplot2::geom_vline(xintercept=0, colour="grey60", linetype='dashed')+
    ggplot2::geom_point(ggplot2::aes(colour=as.numeric(ID_delta)), alpha=0.6)+
    ggplot2::scale_colour_viridis_c("Number of \n observed dimensions",
                                    option="inferno",
                                    guide=ggplot2::guide_colourbar(title.position = "top"))+
    ggplot2::labs(x = expression(paste(Delta, "Sub-optimal species abundance")),
                  y = expression(paste(Delta, " Mean performance")))+
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   axis.text = ggplot2::element_text(size=16),
                   legend.direction="horizontal")
  
  ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_abund_suboptimals_vs_delta_perf_{mortality}_{fecundity}.png")),
                  width=fig_width, height=fig_width, units="cm", dpi=300)
}

# Mean performance vs. sub-optimal species abundance

mortality <- "prop"
fecundity <- "abund"

Perf_vs_suboptimal_abundance <- Comparison_suboptimals[which(Comparison_suboptimals$Mortality==mortality&Comparison_suboptimals$Fecundity==fecundity),]
Perf_vs_suboptimal_abundance <- Perf_vs_suboptimal_abundance[Perf_vs_suboptimal_abundance$Mod!="Perf_know",]

p1 <- ggplot2::ggplot(data=Perf_vs_suboptimal_abundance, ggplot2::aes(x=Abund_suboptimals, y=perf_tot))+
  ggplot2::geom_point(ggplot2::aes(colour=as.factor(Nb_obs)), alpha=0.6)+
  ggplot2::scale_colour_viridis_d("Number of observed dimensions", option="inferno")+
  ggplot2::labs(x = "Abundance of sub-optimal species",
                y = "Mean performance")+
  ggplot2::ggtitle(paste0("Mortality ", mortality))+
  ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size=7))

ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Abund_suboptimals_vs_delta_perf_{mortality}_{fecundity}.png")),
                width=fig_width, height=fig_width, units="cm", dpi=300)

p2 <- ggplot2::ggplot(data=Perf_vs_suboptimal_abundance, ggplot2::aes(x=Abund_suboptimals, y=perf_tot))+
  ggplot2::geom_point(ggplot2::aes(colour=as.factor(Mod)), alpha=0.6)+
  ggplot2::scale_colour_manual("Model", values=c("darkblue","gold"))+
  ggplot2::labs(x = "Abundance of sub-optimal species",
                y = "Mean performance")+
  ggplot2::ggtitle(paste0("Mortality ", mortality))+
  ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size=7))

p <- ggpubr::ggarrange(p1, p2, nrow=1, ncol=2, labels=c("A", "B"))

ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Abund_suboptimals_vs_delta_perf_{mortality}_{fecundity}.png")),
                width=fig_width, height=fig_width, units="cm", dpi=300)

# Repeatability within the same model

nrep <- 10

Percentage_similarity_within <- data.frame(
  Mortality=character(),
  Fecundity=character(),
  Mod=character(),
  Nb_obs=integer(),
  Seed=integer(),
  Seed_r=character(),
  PS=numeric()
)

combi_seed_r <- t(combn(c(1:10), 2))

for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")){
    for (seed in 1:nrep) {
      
      for (k in 1:nrow(combi_seed_r)){
        seed_r_1 <- combi_seed_r[k,1]
        seed_r_2 <- combi_seed_r[k,2]
        load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r_1}_Abundances.RData")))
        abund_1 <- Abundances[ngen,]
        load(here::here("outputs", "outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r_2}_Abundances.RData")))
        abund_2 <- Abundances[ngen,]
        
        Percentage_similarity_within <- rbind(Percentage_similarity_within,
                                              data.frame(
                                                Mortality = mortality,
                                                Fecundity = fecundity,
                                                Mod="Perf_know",
                                                Nb_obs = NA,
                                                Seed = seed,
                                                Seed_r = paste0(seed_r_1, "-", seed_r_2),
                                                PS = percentage_similarity(abund_1, abund_2)
                                              ))
      }
      
      for (mod in c("Part_know", "Part_know_IV")){
        for(nb_obs in 0:n_axes){
          for (k in 1:nrow(combi_seed_r)){
            seed_r_1 <- combi_seed_r[k,1]
            seed_r_2 <- combi_seed_r[k,2]
            load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r_1}_Abundances.RData")))
            abund_1 <- Abundances[ngen,]
            load(here::here("outputs", "outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r_2}_Abundances.RData")))
            abund_2 <- Abundances[ngen,]
            
            Percentage_similarity_within <- rbind(Percentage_similarity_within,
                                                  data.frame(
                                                    Mortality = mortality,
                                                    Fecundity = fecundity,
                                                    Mod=mod,
                                                    Nb_obs = nb_obs,
                                                    Seed = seed,
                                                    Seed_r = paste0(seed_r_1, "-", seed_r_2),
                                                    PS = percentage_similarity(abund_1, abund_2)
                                                  ))
          }
        }
      }
    }
  }
}
save(Percentage_similarity_within, file = here::here("outputs", "Comparison", "Percentage_similarity_within.RData"))

mortality <- "prop"
fecundity <- "abund"

data_figure <- Percentage_similarity_within[which(Percentage_similarity_within$Mortality==mortality&Percentage_similarity_within$Fecundity==fecundity),]
data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- 16
data_figure$Nb_obs <- as.integer(data_figure$Nb_obs)
data_1 <- data_figure[data_figure$Mod!="Perf_know",]
data_2 <- data_figure[data_figure$Mod=="Perf_know",]
data_1_A <- data_1[data_1$Mod!="Part_know",]
data_1_B <- data_1[data_1$Mod!="Part_know_IV",]

color_text <- c(rep("grey20", length(0:n_axes)), "darkred")

pA <- ggplot2::ggplot(data=data_1_A, ggplot2::aes(x=Nb_obs, y=PS))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::scale_colour_viridis_d()+
  ggplot2::geom_boxplot(alpha=0.6, colour="black", ggplot2::aes(group=Nb_obs))+
  ggplot2::geom_jitter(data=data_2, ggplot2::aes(x=Nb_obs, y=PS, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::geom_boxplot(data=data_2, ggplot2::aes(x=Nb_obs, y=PS), alpha=0.6, colour="darkred")+
  ggplot2::scale_x_continuous(breaks=c(0:(n_axes+1)), labels=c(0:n_axes, "PK"))+
  ggplot2::labs(title="IK with uIV",
                x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                y = "Similarity within the model")+
  ggplot2::theme(text = ggplot2::element_text(size = 18),
                 plot.title = ggplot2::element_text(size = 18),
                 axis.title = ggplot2::element_text(size=18),
                 axis.text = ggplot2::element_text(size=13),
                 axis.text.x = ggplot2::element_text(colour = color_text),
                 axis.text.y = ggplot2::element_text(colour = "grey20"),
                 legend.position = "none")

pB <- ggplot2::ggplot(data=data_1_B, ggplot2::aes(x=Nb_obs, y=PS))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::scale_colour_viridis_d()+
  ggplot2::geom_boxplot(alpha=0.6, colour="black", ggplot2::aes(group=Nb_obs))+
  ggplot2::geom_jitter(data=data_2, ggplot2::aes(x=Nb_obs, y=PS, colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::geom_boxplot(data=data_2, ggplot2::aes(x=Nb_obs, y=PS), alpha=0.6, colour="darkred")+
  ggplot2::scale_x_continuous(breaks=c(0:(n_axes+1)), labels=c(0:n_axes, "PK"))+
  ggplot2::labs(title="IK without uIV",
                x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                y = "Similarity within the model")+
  ggplot2::theme(text = ggplot2::element_text(size = 18),
                 plot.title = ggplot2::element_text(size = 18),
                 axis.title = ggplot2::element_text(size=18),
                 axis.text = ggplot2::element_text(size=13),
                 axis.text.x = ggplot2::element_text(colour = color_text),
                 axis.text.y = ggplot2::element_text(colour = "grey20"),
                 legend.position = "none")

p <- ggpubr::ggarrange(pA, pB, ncol=1, nrow=2, labels=c("A", "B"))


data_figure <- Percentage_similarity_within[which(Percentage_similarity_within$Mortality==mortality&Percentage_similarity_within$Fecundity==fecundity),]
data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- 16

p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=PS))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed), shape=as.factor(Mod), group=as.factor(Mod)), alpha=0.6, position=ggplot2::position_jitterdodge(jitter.width=0.5))+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(ggplot2::aes(colour=as.factor(Mod)), alpha=0.6)+
  ggplot2::scale_x_discrete(breaks=c(0:(n_axes+1)), labels=c(0:n_axes, "PK"))+
  ggplot2::scale_colour_manual(values=c("black", "grey50", "darkred"))+
  ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                y = "Similarity within the model")+
  ggplot2::theme(text = ggplot2::element_text(size = 18),
                 plot.title = ggplot2::element_text(size = 18),
                 axis.title = ggplot2::element_text(size=18),
                 axis.text = ggplot2::element_text(size=13),
                 axis.text.x = ggplot2::element_text(colour = color_text),
                 axis.text.y = ggplot2::element_text(colour = "grey20"),
                 legend.position = "none")

ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Percentage_similarity_within_{mortality}_{fecundity}.png")),
                width=fig_width, height=fig_width/1.5, units="cm", dpi=300)

# Results of the Imperfect knowledge model without uIV #

data_figure <- Species_all[which(Species_all$Mortality==mortality&Species_all$Fecundity==fecundity),]
data_figure <- data_figure[data_figure$Mod!="Part_know_IV",]
data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
data_figure$Nb_obs <- factor(data_figure$Nb_obs,
                             levels = c(as.character(c(0:n_axes)), "PK"))

data_figure$Boxplot_colour <- rep(0, nrow(data_figure))
data_figure[which(data_figure$Mod=="Perf_know"),]$Boxplot_colour <- 1
color_text <- c(rep("grey20", length(0:n_axes)), "darkred")

p1 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=N_sp))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
  ggplot2::scale_colour_manual(values=c("black", "darkred"))+
  ggplot2::labs(x = "Number of observed dimensions",
                y = "Species richness without uIV")+
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 axis.text = ggplot2::element_text(size=14),
                 axis.text.x = ggplot2::element_text(colour = color_text),
                 axis.text.y = ggplot2::element_text(colour = "grey20"),
                 legend.position = "none")

p2 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Hill_Shannon))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
  ggplot2::scale_colour_manual(values=c("black", "darkred"))+
  ggplot2::labs(x = "Number of observed dimensions",
                y = "Hill-Shannon index without uIV")+
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 axis.text = ggplot2::element_text(size=14),
                 axis.text.x = ggplot2::element_text(colour = color_text),
                 axis.text.y = ggplot2::element_text(colour = "grey20"),
                 legend.position = "none")

plot_level_explanation <- ggplot2::ggplot(data_level_explanation, ggplot2::aes(Nb_obs, Mean_explanation))+
  ggplot2::geom_blank()+
  ggplot2::theme_classic()+
  ggplot2::scale_x_continuous(name="Level of explained variance", breaks=data_level_explanation$Nb_obs, labels=labels)+
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 axis.text = ggplot2::element_text(size=13),
                 axis.line.x = ggplot2::element_line(arrow = grid::arrow(length=grid::unit(0.15, "inches")), colour="deeppink3"),
                 axis.ticks.x=ggplot2::element_line(colour="deeppink3"),
                 axis.text.x=ggplot2::element_text(colour="deeppink3"),
                 axis.line.y=ggplot2::element_blank(),
                 axis.text.y=ggplot2::element_blank(),
                 axis.ticks.y=ggplot2::element_blank(),
                 axis.title.y=ggplot2::element_blank(),
                 panel.grid.minor.y=ggplot2::element_blank(),
                 panel.grid.major.y=ggplot2::element_blank(),
                 plot.margin = ggplot2::margin(l = 36,
                                               r = 25,
                                               t = 0,
                                               b = 0.5))

p <- ggpubr::ggarrange(p1, p2,
                       nrow=1, ncol=2, align="v", labels=c("A", "B"))

arrange_arrow <- ggpubr::ggarrange(plot_level_explanation, plot_level_explanation,
                                   nrow=1, ncol=2, align="v")

p_S_SR_Hill_Shannon <- ggpubr::ggarrange(p, arrange_arrow,
                                    nrow=2, ncol=1, heights=c(10,1))

ggplot2::ggsave(p_S_SR_Hill_Shannon, filename=here::here("outputs", "Comparison", glue::glue("Results_1_IK_without_uIV_{mortality}_{fecundity}.png")),
                width=fig_width, height=fig_width/2, units="cm", dpi=300)


data_figure <- Percentage_similarity[which(Percentage_similarity$Mortality==mortality&Percentage_similarity$Fecundity==fecundity),]
data_figure <- data_figure[data_figure$Mod!="Part_know_IV",]
data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
data_figure$Nb_obs <- factor(data_figure$Nb_obs,
                             levels = c(as.character(c(0:n_axes)), "PK"))

data_figure$Boxplot_colour <- rep(0, nrow(data_figure))
data_figure[which(data_figure$Mod=="Perf_know"),]$Boxplot_colour <- 1
color_text <- c(rep("grey20", length(0:n_axes)), "darkred")

p1 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=PS))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
  ggplot2::scale_colour_manual(values=c("black", "darkred"))+
  ggplot2::labs(x = "Number of observed dimensions",
                y = "Similarity between PK and IK without uIV")+
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 axis.text = ggplot2::element_text(size=14),
                 axis.text.x = ggplot2::element_text(colour = color_text),
                 axis.text.y = ggplot2::element_text(colour = "grey20"),
                 legend.position = "none")

data_figure <- Perf_on_sites[which(Perf_on_sites$Mortality==mortality&Perf_on_sites$Fecundity==fecundity),]
data_figure <- data_figure[data_figure$Mod!="Part_know_IV",]
data_figure[which(data_figure$Mod=="Perf_know"),]$Nb_obs <- "PK"
data_figure$Nb_obs <- factor(data_figure$Nb_obs,
                             levels = c(as.character(c(0:n_axes)), "PK"))

data_figure$Boxplot_colour <- rep(0, nrow(data_figure))
data_figure[which(data_figure$Mod=="Perf_know"),]$Boxplot_colour <- 1
color_text <- c(rep("grey20", length(0:n_axes)), "darkred")

p2 <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(Nb_obs), y=Perf))+
  ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
  ggplot2::scale_colour_viridis_d()+
  ggnewscale::new_scale("colour")+
  ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(colour=as.factor(Boxplot_colour)))+
  ggplot2::scale_colour_manual(values=c("black", "darkred"))+
  ggplot2::labs(x = "Number of observed dimensions",
                y = "Mean performance without uIV")+
  ggplot2::theme(text = ggplot2::element_text(size = 16),
                 axis.text = ggplot2::element_text(size=14),
                 axis.text.x = ggplot2::element_text(colour = color_text),
                 axis.text.y = ggplot2::element_text(colour = "grey20"),
                 legend.position = "none")

p <- ggpubr::ggarrange(p1, p2,
                       nrow=1, ncol=2, align="v", labels=c("A", "B"))

arrange_arrow <- ggpubr::ggarrange(plot_level_explanation, plot_level_explanation,
                                   nrow=1, ncol=2, align="v")

p_S_PS_Perf <- ggpubr::ggarrange(p, arrange_arrow,
                                 nrow=2, ncol=1, heights=c(10,1))

ggplot2::ggsave(p_S_PS_Perf, filename=here::here("outputs", "Comparison", glue::glue("Results_2_IK_without_uIV_{mortality}_{fecundity}.png")),
                width=fig_width, height=fig_width/2, units="cm", dpi=300)
