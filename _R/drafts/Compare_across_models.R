perf_know <- TRUE
IV <- FALSE

source(file = here::here("_R", "Basic_parameters.R"))

#Dependence between the number of seeds and the abundance of species
nb_seeds_dep_abund<-TRUE

#Probability of mortality as a function of species performance
mortality_fixed<-TRUE

compare_models(nb_obs_axes, Seeds, nrep, nsp, ngen, nsite_side, n_axes)

Compare_IV_axis_nb(Seeds, nsp, nb_obs_axes)

#Dependence between the number of seeds and the abundance of species
nb_seeds_dep_abund<-TRUE

#Probability of mortality as a function of species performance
mortality_fixed<-FALSE

source(file = here::here("_R", "Basic_parameters.R"))

compare_models(nb_obs_axes, Seeds, nrep, nsp, ngen, nsite_side, n_axes)

Compare_IV_axis_nb(Seeds, nsp, nb_obs_axes)

#Dependence between the number of seeds and the abundance of species
nb_seeds_dep_abund<-FALSE

#Probability of mortality as a function of species performance
mortality_fixed<-FALSE

source(file = here::here("_R", "Basic_parameters.R"))

compare_models(nb_obs_axes, Seeds, nrep, nsp, ngen, nsite_side, n_axes)

Compare_IV_axis_nb(Seeds, nsp, nb_obs_axes)

#Dependence between the number of seeds and the abundance of species
nb_seeds_dep_abund<-FALSE

#Probability of mortality as a function of species performance
mortality_fixed<-TRUE

source(file = here::here("_R", "Basic_parameters.R"))

compare_models(nb_obs_axes, Seeds, nrep, nsp, ngen, nsite_side, n_axes)

Compare_IV_axis_nb(Seeds, nsp, nb_obs_axes)

