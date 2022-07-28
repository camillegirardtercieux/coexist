# Compute multidimensional semivariance #

#install archive of RandomFields
#install.packages("geoR", contriburl = "http://www.leg.ufpr.br/geoR")

Correlation_env_sp <- data.frame(
  Mortality = numeric(),
  Fecundity = numeric(),
  Seed = numeric(),
  Seed_r = numeric(),
  Mod = character(),
  Nb_obs = numeric(),
  Correlation = numeric()
  )

for (simu in c(1:nrow(Simulations))){
  print(paste("Simu", simu))

  mortality <- Simulations[simu,1]
  fecundity <- Simulations[simu,2]
  seed <- Simulations[simu,3]
  seed_r <- Simulations[simu,4]

  #Load the environment and species optima once per configuration (seed)
  if(seed_r==1){
    load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_sites.RData")))
    load(here::here("outputs_cluster", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_niche_optimum.RData")))
  }

  for (mod in c("Perf_know", "Part_know", "Part_know_IV")){

    if(mod == "Perf_know"){

      load(here::here("outputs_cluster", glue::glue("{mod}_0_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
      sp_XY <- data.frame(raster::rasterToPoints(raster::raster(community_end)))
      names(sp_XY) <- c("x", "y", "sp")
      vario_sp <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)

      semivar_multidim <- compute_semivar_multidim(sites, n_axes, niche_optimum, sp_XY, vario_sp, nsp, community_end)
      semivar_multidim$Vario_sp_geoR <- vario_sp$u
      semivar_multidim$Distance <- vario_sp$bins.lim[-length(vario_sp$bins.lim)]
      semivar_multidim$Sample_size <- vario_sp$n

      semivar_multidim<-semivar_multidim%>%
        filter(Sample_size>500)

      m <- lm(semivar_multidim$Vario_sp ~ semivar_multidim$Vario_env)

      Correlation_env_sp_temp <- data.frame(
        Mortality = mortality,
        Fecundity = fecundity,
        Seed = seed,
        Seed_r = seed_r,
        Mod = mod,
        Nb_obs = NA,
        Correlation = round(sqrt(summary(m)$r.squared), digits = 2)
      )

      Correlation_env_sp <- rbind(Correlation_env_sp, Correlation_env_sp_temp)

    }else{

      for (nb_obs in c(0:15)){

        load(here::here("outputs_cluster", glue::glue("{mod}_{nb_obs}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
        sp_XY <- data.frame(raster::rasterToPoints(raster::raster(community_end)))
        names(sp_XY) <- c("x", "y", "sp")
        vario_sp <- geoR::variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)

        semivar_multidim <- compute_semivar_multidim(sites, n_axes, niche_optimum, sp_XY, vario_sp, nsp, community_end)
        semivar_multidim$Vario_sp_geoR <- vario_sp$u
        semivar_multidim$Distance <- vario_sp$bins.lim[-length(vario_sp$bins.lim)]
        semivar_multidim$Sample_size <- vario_sp$n

        semivar_multidim<-semivar_multidim%>%
          filter(Sample_size>500)

        m <- lm(semivar_multidim$Vario_sp ~ semivar_multidim$Vario_env)

        Correlation_env_sp_temp <- data.frame(
          Mortality = mortality,
          Fecundity = fecundity,
          Seed = seed,
          Seed_r = seed_r,
          Mod = mod,
          Nb_obs = nb_obs,
          Correlation = round(sqrt(summary(m)$r.squared), digits = 2)
        )

        Correlation_env_sp <- rbind(Correlation_env_sp, Correlation_env_sp_temp)

      }#for n_obs
    }#else Part_know
  }#for mod
}#for simu

save(Correlation_env_sp, file=here::here("outputs", "Comparison", "Correlation_env_sp.RData"))

load(file=here::here("outputs", "Comparison", "Correlation_env_sp.RData"))

#One plot per Mortality * Fecundity option
for (mortality in c("fixed", "prop", "stocha", "stocha_basal")){
  for (fecundity in c("abund", "fixed")){

    data_figure <- Correlation_env_sp[which(Correlation_env_sp$Mortality==mortality&Correlation_env_sp$Fecundity==fecundity),]

    Summary_correlation_env_sp<-data_figure%>%
      dplyr::group_by(Mod, Nb_obs)%>%
      dplyr::mutate(Mean_corr=mean(Correlation, na.rm=TRUE), Sd=sd(Correlation, na.rm=TRUE))%>%
      dplyr::slice(1)%>%
      dplyr::ungroup()%>%
      dplyr::select(-Correlation, -Seed, -Seed_r)

    save(Summary_correlation_env_sp, file=here::here("outputs", "Comparison", glue::glue("Mean_correlation_env_sp_{mortality}_{fecundity}.RData")))

    data_figure <- data_figure[data_figure$Mod!="Perf_know",]

    #Compute delta between with and without
    data_figure <- data_figure%>%
      dplyr::group_by(Seed, Seed_r, Nb_obs)%>%
      dplyr::mutate(Delta = Correlation - dplyr::lag(Correlation) , ID_delta=Nb_obs)%>%
      dplyr::filter(is.na(Delta)==FALSE)%>%
      dplyr::ungroup()%>%
      dplyr::select(Seed, Seed_r, ID_delta, Delta)

    p <- ggplot2::ggplot(data=data_figure, ggplot2::aes(x=as.factor(ID_delta), y=Delta))+
      ggplot2::geom_hline(yintercept=0, colour="grey60", linetype='dashed')+
      ggplot2::geom_jitter(ggplot2::aes(colour=as.factor(Seed)), alpha=0.6, height=0, width=0.3, shape=16)+
      ggplot2::geom_boxplot(alpha=0.6, ggplot2::aes(group=ID_delta))+
      ggplot2::scale_colour_viridis_d()+
      ggplot2::labs(x = expression(paste("Number of observed dimensions ( ~ ", frac(sIV,uIV), " )")),
                    y = expression(paste(Delta, " environment-species correlation")))+
      ggplot2::theme(axis.title = ggplot2::element_text(size = 10),
                     axis.text = ggplot2::element_text(size=7),
                     legend.position = "none")

    ggplot2::ggsave(p, filename=here::here("outputs", "Comparison", glue::glue("Delta_corr_env_sp_{mortality}_{fecundity}.png")),
                    width=fig_width, height=fig_width/2, units="cm", dpi=300)
  }
}
