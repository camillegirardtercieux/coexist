launch_model <- function(mortality, fecundity, seed, seed_r){
  
  if (mortality =="fixed"){
    mortality_fixed<-TRUE
    mortality_proportion<-FALSE
    mortality_stocha<-FALSE
    mortality_stocha_basal<-FALSE
  }
  
  if (mortality =="prop"){
    mortality_fixed<-FALSE
    mortality_proportion<-TRUE
    mortality_stocha<-FALSE
    mortality_stocha_basal<-FALSE
  }
  
  if (mortality =="stocha"){
    mortality_fixed<-FALSE
    mortality_proportion<-FALSE
    mortality_stocha<-TRUE
    mortality_stocha_basal<-FALSE
  }
  
  if (mortality =="stocha_basal"){
    mortality_fixed<-FALSE
    mortality_proportion<-FALSE
    mortality_stocha<-TRUE
    mortality_stocha_basal<-TRUE
  }
  
  if(fecundity=="fixed"){
    nb_seeds_dep_abund<-FALSE
    nb_seeds_indep_abund<-TRUE
  }
  
  if(fecundity=="abund"){
    nb_seeds_dep_abund<-TRUE
    nb_seeds_indep_abund<-FALSE
  }
  
  source(file = here::here("_R", "basic_parameters.R"), local=TRUE)
  
  
  # =========================
  # Landscape and environment
  # =========================
  
  if(perf_know==TRUE){
    generate_environment(n_axes=n_axes, nsite_side=nsite_side, seed_env=Seeds[seed], mod=mod, n_observed_axes=n_observed_axes, mortality=mortality, fecundity=fecundity, seed=seed, seed_r=seed_r)
    load(paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_sites.RData")))
    load(paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_env.RData")))
  } else {
    if(n_observed_axes>0){
      Obs_env <- sites[,1:n_observed_axes]
    }
  }
  
  # =========================================
  # Species optima
  # =========================================
  
  if(perf_know==TRUE){
    niche_optimum <- NULL
    niche_optimum <- generate_species_optima(randomOptSp=randomOptSp, niche_width=niche_width, nsp=nsp, env=env, n_axes=n_axes, seed_sp=Seeds[seed], mod=mod, n_observed_axes=n_observed_axes, mortality=mortality, fecundity=fecundity, seed=seed, seed_r=seed_r)
    load(paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_niche_optimum.RData")))
  } else {
    
    Inferred_species_parameters_mat <-list()
    
    for(k in 1:ncol(Inferred_species_parameters)){
      Inferred_species_parameters_mat[[k]] <- matrix(rep(Inferred_species_parameters[,k],each=nsite), ncol=nsp)
    }
    
    if(n_observed_axes>0){
      
      Obs_env_mat <- list()
      
      if(is.null(ncol(Obs_env))==FALSE){
        
        for(k in 1:ncol(Obs_env)){
          Obs_env_mat[[k]] <- matrix(rep(Obs_env[,k], nsp), ncol=nsp)
          Obs_env_mat[[k+ncol(Obs_env)]] <- matrix(rep((Obs_env[,k])^2, nsp), ncol=nsp)
        }
      }
      
      else{
        Obs_env_mat[[1]] <- matrix(rep(Obs_env, nsp), ncol=nsp)
        Obs_env_mat[[2]] <- matrix(rep(Obs_env^2, nsp), ncol=nsp)
      }
    }
  }
  
  # =========================================
  # Species performance
  # =========================================
  
  # Matrix of species performance on each site (distances)
  # Sites in rows, Species in columns
  if(perf_know==TRUE){
    dist_E_Sp <- dist_Site_Sp(as.matrix(sites), as.matrix(niche_optimum))
    dprim_E_Sp <- (dist_E_Sp-mean(dist_E_Sp))/sd(dist_E_Sp)
    perf_E_Sp <- -dprim_E_Sp
    save(perf_E_Sp, file=paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_perf_E_Sp.RData")))
  } else {
    
    perf_Sp_mean <- Inferred_species_parameters_mat[[1]]
    
    if(n_observed_axes>0){
      
      for(k in 1:length(Obs_env_mat)){
        perf_Sp_mean <- perf_Sp_mean + Inferred_species_parameters_mat[[k+1]]*Obs_env_mat[[k]]
      }
    }
    save(perf_Sp_mean, file = paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_perf_Sp_mean.RData")))
  }
  
  # =========================================
  # Species mortality probability
  # =========================================
  
  if(perf_know==TRUE){
    if(mortality_stocha==TRUE){
      # Probability of dying of each species on each site
      if(mortality_stocha_basal==TRUE){
        mortality_E_Sp <- inv_logit(logit(theta) + b * perf_E_Sp)
      }
      if(mortality_stocha_basal==FALSE){
        mortality_E_Sp <- inv_logit(b*perf_E_Sp)
      }
      
      save(mortality_E_Sp, file = paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_mortality_E_Sp.RData")))
      
    }#end condition on stochastic mortality
    
  }#end condition on perfect knowledge
  
  if(perf_know==FALSE){
    
    if(mortality_stocha==TRUE){
      if(mortality_stocha_basal==TRUE){
        mortality_Sp_mean <- inv_logit(logit(theta) + b * (perf_Sp_mean))
      }
      if(mortality_stocha_basal==FALSE){
        mortality_Sp_mean <- inv_logit(b*perf_Sp_mean)
      }
    }#end condition on stochastic mortality
  }#end condition on partial knowledge
  
  # =========================================
  # Repetitions
  # =========================================
  
  #Abundance matrices
  Abundances <- matrix(NA, ncol=nsp, nrow=ngen)
  # Mean mortality rate in the community
  theta_comm <- c()
  mortality_rate <- c()
  
  abund <- matrix(NA, ncol=nsp, nrow=ngen)
  #used if nb_seeds_dep_abund==TRUE 
  abund_after_mortality <- NULL
  
  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  if(start_full_landscape==TRUE){
    # Draw species at random in the landscape (one individual per site)
    set.seed(Seeds[seed_r])
    sp <- sample(1:nsp, size=nsite, replace=TRUE)
    community <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
  }
  if(start_one_ind_per_species==TRUE){
    # One individual per species distributed randomly over the grid
    set.seed(Seeds[seed_r])
    sites_start <- sample(1:nsite, size=nsp, replace=FALSE)
    set.seed(Seeds[seed_r])
    sp_start<- sample(1:nsp, size=nsp, replace=FALSE)
    community <- matrix(0, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    community_rast <- raster::raster(community)
    community_rast[sites_start] <- sp_start
    community <- raster::as.matrix(community_rast)
  }
  if(start_ten_ind_per_species==TRUE){
    # Ten individuals per species distributed randomly over the grid
    set.seed(Seeds[seed_r])
    sites_start <- sample(1:nsite, size=nsp*10, replace=FALSE)
    sp_start<- rep(1:nsp, each=10)
    community <- matrix(0, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    community_rast <- raster::raster(community)
    community_rast[sites_start] <- sp_start
    community <- raster::as.matrix(community_rast)
  }
  
  # Plot the community at the start of the first repetition
  if (seed_r==1) {
    save(community, file=paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_community_start.RData")))
  }
  
  if(perf_know==FALSE&&IV==TRUE){
    # Adding endogenous IV to compute mortality for the first generation
    epsilon <- rnorm(nsite*nsp, 0, sd=sqrt(rep(V_intra$V, each=nsite)))
    perf_ind <- perf_Sp_mean + epsilon
  }
  
  # -----------------------------------------
  # Dynamics
  # -----------------------------------------
  
  # Simulating generation
  for (g in 1:ngen) {
    
    # Species richness
    abund[g, ] <- table(factor(as.vector(community), levels=1:nsp))
    
    # ******************
    # Mortality
    # ******************
    
    #Total abundance in the community, used for mortality rate as a proportion of abundance
    abund_tot <- length(community[community!=0])
    
    if(mortality_stocha==TRUE){
      
      if(perf_know==FALSE&&IV==TRUE){
        
        if(mortality_stocha_basal==TRUE){
          # Probability of dying of each species on each site
          mortality_ind <- inv_logit(logit(theta) + b * (perf_ind))
          
          # Mean mortality correction
          epsilon_mat <- matrix(epsilon, ncol=nsp)
          theta_var <- inv_logit(logit(theta) + b * epsilon_mat)
          diff <- mean(theta_var)-theta
          mortality_ind <- mortality_ind - diff
          # /!\ can modify mean!
          mortality_ind[mortality_ind<0] <- 0
        }#end condition on basal mortality
        
        if(mortality_stocha_basal==FALSE){
          mortality_ind <- inv_logit(b*perf_ind)
        }
        
      } # end condition on IV
      
      # Mortality rate on each site
      # w0 for vacant sites
      theta_site <- rep(NA, nsite)
      w0 <- (as.vector(t(community))==0)
      theta_site[w0] <- 0
      
      if(perf_know==TRUE&&IV==FALSE){
        if(length(c(theta_site[!w0]))>1){
          theta_site[!w0] <- diag(mortality_E_Sp[!w0, as.vector(t(community))[!w0]])
        }
        if(length(c(theta_site[!w0]))==1){
          theta_site[!w0] <- mortality_E_Sp[!w0, as.vector(t(community))[!w0]]
        }
      }# end condition perf_know TRUE
      
      if(perf_know==FALSE&&IV==FALSE){
        if(length(c(theta_site[!w0]))>1){
          theta_site[!w0] <- diag(mortality_Sp_mean[!w0, as.vector(t(community))[!w0]])
        }
        if(length(c(theta_site[!w0]))==1){
          theta_site[!w0] <- mortality_Sp_mean[!w0, as.vector(t(community))[!w0]]
        }
      }# end condition IV FALSE
      
      if(perf_know==FALSE&&IV==TRUE){
        if(length(c(theta_site[!w0]))>1){
          theta_site[!w0] <- diag(mortality_ind[!w0, as.vector(t(community))[!w0]])
        }
        if(length(c(theta_site[!w0]))==1){
          theta_site[!w0] <- mortality_ind[!w0, as.vector(t(community))[!w0]]
        }
      }# end condition IV TRUE
      
      # Mortality events
      if(mortality_stocha_basal==TRUE){
        mort_ev <- rbinom(nsite, size=1, prob=theta_site)
      }
      if(mortality_stocha_basal==FALSE){
        ind_dead <- sample(x=1:nsite, size=round(0.01*abund_tot), replace=FALSE, prob=theta_site)
        mort_ev <- rep(0, nsite)
        mort_ev[ind_dead] <- 1
      }
      mortality_matrix <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
      
    }#end condition on mortality stochastic
    
    if(mortality_stocha==FALSE){
      # No stochasticity: the 10 less performing individuals of the community OR a fixed proportion of the total abundance die each generation
      
      perf_present <- rep(NA, nsite)
      w0 <- (as.vector(t(community))==0)
      perf_present[w0] <- NA
      
      if(perf_know==TRUE&&IV==FALSE){
        #performance of the present individuals
        #empty sites are not taken into account
        perf_present[!w0] <- diag(perf_E_Sp[!w0, as.vector(t(community))[!w0]])
      }
      
      if(perf_know==FALSE&&IV==FALSE){
        #performance of the present individuals
        perf_present[!w0] <- diag(perf_Sp_mean[!w0, as.vector(t(community))[!w0]])
      }
      
      if(perf_know==FALSE&&IV==TRUE){
        #performance of the present individuals
        perf_present[!w0] <- diag(perf_ind[!w0, as.vector(t(community))[!w0]])
      }
      
      mort_ev <- rep(0, length(perf_present))
      #identify the 10 less performant present individuals and kill them
      if(mortality_fixed==TRUE){
        if(length(which(perf_present<=sort(perf_present)[10]))<=10){
          mort_ev[which(perf_present<=sort(perf_present)[10])] <- 1 #sort ignores NA
        }else{
          keep_sp <- which(perf_present<sort(perf_present)[10])
          sample_sp <- sample(which(perf_present==sort(perf_present)[10]), 10-length(keep_sp))
          mort_ev[c(keep_sp, sample_sp)] <- 1
        }
      }
      #identify the 0.01*total abundance less performing present individuals and kill them
      if(mortality_proportion==TRUE){
        if(length(which(perf_present<=sort(perf_present)[round(abund_tot*0.01)]))<=round(abund_tot*0.01)){
          mort_ev[which(perf_present<=sort(perf_present)[round(abund_tot*0.01)])] <- 1 #sort ignores NA
        }else{
          keep_sp <- which(perf_present<sort(perf_present)[round(abund_tot*0.01)])
          sample_sp <- sample(which(perf_present==sort(perf_present)[round(abund_tot*0.01)]), round(abund_tot*0.01)-length(keep_sp))
          mort_ev[c(keep_sp, sample_sp)] <- 1
        }
      }
      mortality_matrix <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    }
    
    # Number of deaths
    n_mort <- sum(mort_ev)
    mortality_rate[g] <- n_mort
    
    if(n_mort!=0){
      
      # Update community
      community[mortality_matrix==1] <- 0
      
      
      # *********************
      # Fecundity/Recruitment
      # *********************
      
      # Species present in the community
      sp_present <- sort(unique(community[community!=0]))
      nsp_present <- length(sp_present)
      
      # Vacant sites
      community_rast <- raster::raster(community)
      sites_vacant <- which(raster::values(community_rast)==0)
      nsite_vacant <- length(sites_vacant)
      
      if(perf_know==FALSE&&IV==TRUE){
        # New individual effects for potential new individuals (one per site and species)
        epsilon <- rnorm(nsite*nsp, 0, sd=sqrt(rep(V_intra$V, each=nsite)))
        perf_ind_pot <- perf_Sp_mean + epsilon
      }
      
      if(nb_seeds_dep_abund==FALSE){
        
        if(perf_know==TRUE&&IV==FALSE){
          
          # Performance of species on vacant sites
          dist_E_Sp_vacant <- dist_E_Sp[sites_vacant, ]
          
          # Identify the present species with the highest performance on vacant sites
          if(!is.null(nrow(dist_E_Sp_vacant))){
            new_ind <- apply(dist_E_Sp_vacant, 1, high_perf_sp, sp_pres=sp_present)
          }else{new_ind <- high_perf_sp(dist=dist_E_Sp_vacant, sp_pres=sp_present)}
          
          # Recruitment
          community_rast[sites_vacant] <- new_ind
          
        } # end condition perfect knowledge and no IV
        
        if(perf_know==FALSE&&IV==FALSE){
          
          # Identify the present species with the highest performance on vacant sites
          if(!is.null(nrow(-matrix(perf_Sp_mean[sites_vacant, ], ncol=nsp)))){
            sp_high_perf <- apply(-matrix(perf_Sp_mean[sites_vacant, ], ncol=nsp), 1, high_perf_sp, sp_pres=sp_present)
          }else{sp_high_perf <- high_perf_sp(dist=-matrix(perf_Sp_mean[sites_vacant, ], ncol=nsp), sp_pres=sp_present)}
          
          
          # Recruitment
          community_rast[sites_vacant] <- sp_high_perf
          
        } # end condition partial knowledge and no IV
        
        if(perf_know==FALSE&&IV==TRUE){
          
          # Identify the present species with the highest performance on vacant sites (maximum in each line)
          if(!is.null(nrow(-matrix(perf_ind_pot[sites_vacant, ], ncol=nsp)))){
            sp_high_perf <- apply(-matrix(perf_ind_pot[sites_vacant, ], ncol=nsp), 1, high_perf_sp, sp_pres=sp_present)
          }else{sp_high_perf <- high_perf_sp(dist=-matrix(perf_ind_pot[sites_vacant, ], ncol=nsp), sp_pres=sp_present)}
          
          
          # Recruitment
          community_rast[sites_vacant] <- sp_high_perf
          
          # Update performance matrix
          for(k in 1:nsite_vacant){
            perf_ind[sites_vacant[k], sp_high_perf[k]] <- perf_ind_pot[sites_vacant[k], sp_high_perf[k]]
          }
          
        } # end condition partial knowledge + IV
        
      } # end condition no abundance dependence
      
      if (nb_seeds_dep_abund==TRUE){
        
        abund_after_mortality <- as.data.frame(table(factor(as.vector(community), levels=1:nsp)))$Freq
        nb_seeds_sp <- round(abund_after_mortality*fec)
        nb_seeds_tot <- sum(nb_seeds_sp)
        
        # Each seed is dispersed to a random vacant site; several seeds can land on the same site.
        
        sites_of_seeds <- sample(sites_vacant, nb_seeds_tot, replace=TRUE)
        
        # Performance on vacant sites of species which have seeds
        seeds <- data.frame(Species=rep(1:nsp, times=nb_seeds_sp), Sites=sites_of_seeds)
        
        if (nb_seeds_tot>1) {
          if(perf_know==TRUE&&IV==FALSE){
            seeds$Perf <- diag(perf_E_Sp[seeds$Sites, seeds$Species])
          }
          if(perf_know==FALSE&&IV==FALSE){
            seeds$Perf <- diag(perf_Sp_mean[seeds$Sites, seeds$Species])
          }
          if(perf_know==FALSE&&IV==TRUE){
            seeds$Perf <- diag(perf_ind_pot[seeds$Sites, seeds$Species])
          }
        } else {
          if(perf_know==TRUE&&IV==FALSE){
            seeds$Perf <- perf_E_Sp[seeds$Sites, seeds$Species]
          }
          if(perf_know==FALSE&&IV==FALSE){
            seeds$Perf <- perf_Sp_mean[seeds$Sites, seeds$Species]
          }
          if(perf_know==FALSE&&IV==TRUE){
            seeds$Perf <- perf_ind_pot[seeds$Sites, seeds$Species]
          }
        } # end else (one seed)
        
        new_ind <- plyr::ddply(seeds, c("Sites"), f_sp)
        colnames(new_ind) <- c("Sites", "Species")
        
        # Recruitment
        community_rast[new_ind$Sites] <- new_ind$Species
        
        # Update performance matrix with recruited seeds
        if(perf_know==FALSE&&IV==TRUE){
          for(k in 1:nrow(new_ind)){
            perf_ind[new_ind$Sites[k], new_ind$Species[k]] <- perf_ind_pot[new_ind$Sites[k], new_ind$Species[k]]
          }
        }
        
      } # end condition on abundance dependence of number of seeds
      
    }  # end condition on n_mort (contains all cases)
    
    # update community
    community <- raster::as.matrix(community_rast)
    
    # Mean mortality proba in the community
    if(mortality_stocha==TRUE){
      theta_comm[g] <- mean(theta_site[!w0])
    }
    
  } # End ngen
  
  # store final community
  community_end <- community
  Abundances <- abund
  
  save(community_end, file=paste0(directory_writing,"/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_community_end.RData")))
  save(Abundances, file=paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_Abundances.RData")))
  save(mortality_rate, file=paste0(directory_writing ,"/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_mortality_rate.RData")))
  if(mortality_stocha==TRUE){
    save(theta_comm, file=paste0(directory_writing, "/outputs/", glue::glue("{mod}_{n_observed_axes}_{mortality}_{fecundity}_{seed}_{seed_r}_theta_comm.RData")))
  }
}