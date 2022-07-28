## ==============================================================================
## authors      :Ghislain Vieilledent & Camille Girard-Tercieux
## emails       :ghislain.vieilledent@cirad.fr & camillegirardtercieux@gmail.com
## license      :GPLv3
## ==============================================================================

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
s <- as.integer(slurm_arrayid)

directory_reading <- "/home/girardtercieuxc/Chap_2"
directory_writing <- "/lustre/girardtercieuxc"

# To run locally:
# directory_writing <- here::here()
# directory_reading <- here::here()
# add a loop "for (s in 1:nrow(Simulations)){}

setwd(directory_reading)

source(file = paste0(directory_reading, "/_R/call_libraries.R"), local=TRUE)

source(file=here::here("_R", "math_functions.R"), local=TRUE)

source(file=here::here("_R", "generate_environment.R"), local=TRUE)

source(file=here::here("_R", "species_parameters.R"), local=TRUE)

source(file = here::here("_R", "launch_model.R"), local=TRUE)

source(here::here("_R", "infer_IV.R"), local=TRUE)

# Number of observed axes in partial models
nb_obs_axes <- c(0:15)

# Seeds for reproducibility: it controls the environment × species parameters configuration.
# Seeds <- sample(1:10^6, 10)
# save(Seeds, file=here::here("outputs", "Seeds.RData"))

load(here::here("Seeds.RData"))

#source(here::here("_R", "create_simulation_array.R"), local=TRUE)

load(here::here("Array_simulations.RData"))

# ========================================
# Launch perfect knowledge model
# ========================================

perf_know <- TRUE
IV <- FALSE

n_observed_axes <- 0

launch_model(mortality=Simulations[s, 1], fecundity=Simulations[s, 2], seed=as.numeric(Simulations[s, 3]), seed_r=as.numeric(Simulations[s, 4]))

# ========================================
# Infer observed intraspecific variability
# ========================================

for(n_observed_axes in nb_obs_axes){
  infer_IV(n_observed_axes=n_observed_axes, mortality=Simulations[s, 1], fecundity=Simulations[s, 2], seed=as.numeric(Simulations[s, 3]), seed_r=as.numeric(Simulations[s, 4]))
}

# ========================================
# Launch partial knowledge models
# ========================================

perf_know <- FALSE

for(n_observed_axes in nb_obs_axes){
  
  # Without IV
  IV <- FALSE
  
  launch_model(mortality=Simulations[s, 1], fecundity=Simulations[s, 2], seed=as.numeric(Simulations[s, 3]), seed_r=as.numeric(Simulations[s, 4]))
  
  #With IV
  IV <- TRUE
  
  launch_model(mortality=Simulations[s, 1], fecundity=Simulations[s, 2], seed=as.numeric(Simulations[s, 3]), seed_r=as.numeric(Simulations[s, 4]))
}

# # =========================
# # End of file
# # =========================