infer_IV <- function(n_observed_axes, mortality, fecundity, seed, seed_r){
  
  load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_perf_E_Sp.RData")))
  load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_env.RData")))
  load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_sites.RData")))
  load(paste0(directory_writing, "/outputs/", glue::glue("Perf_know_0_{mortality}_{fecundity}_{seed}_{seed_r}_niche_optimum.RData")))
  
  nsp <- ncol(perf_E_Sp)
  n_axes <- length(env)
  
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
  
  # Observed intraspecific variability
  
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
  
  Inferred_species_parameters <- data.frame(matrix(nrow=nsp, ncol=1+2*n_observed_axes))
  Inferred_species_parameters[,1] <- as.vector(lm_fit$coefficients[1:nsp])
  colnames(Inferred_species_parameters)[1]<-"beta_0"
  
  if(n_observed_axes>0){
    
    for(k in 1:(2*n_observed_axes)){
      colnames(Inferred_species_parameters)[k+1] <- glue::glue("beta_{k}")
    }
    
    for(k in 1:n_observed_axes){
      Inferred_species_parameters[,k+1] <- as.vector(c(lm_fit$coefficients[(k*nsp+1):((k+1)*nsp)]))
      Inferred_species_parameters[,k+1+n_observed_axes] <- as.vector(c(lm_fit$coefficients[(nsp*(n_observed_axes+k)+1):(nsp*(n_observed_axes+k+1))]))
    }
  }
  
  save(Inferred_species_parameters, file=paste0(directory_writing, "/outputs/", glue::glue("Inferred_species_parameters_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_{seed_r}.RData")))
  
  load(file=paste0(directory_writing, "/outputs/", glue::glue("Inferred_species_parameters_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_{seed_r}.RData")))
  
  V_intra <- df_perf %>%
    mutate(res=lm_fit$residuals) %>%
    group_by(Species) %>%
    summarise(V=var(res))
  save(V_intra, file = paste0(directory_writing, "/outputs/", glue::glue("V_intra_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_{seed_r}.RData")))
  
  load(paste0(directory_writing, "/outputs/", glue::glue("V_intra_{n_observed_axes}_obs_axes_{mortality}_{fecundity}_{seed}_{seed_r}.RData")))
}